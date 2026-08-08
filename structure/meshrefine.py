"""Utilities for iterative k-mesh refinement based on energy-space error estimates.
   2026 Juntian Peng"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Tuple

import numpy as np

try:  # mpmath provides high-precision polygamma
    import mpmath as mp
except ImportError:  # pragma: no cover - dependency should exist in runtime
    mp = None

from scipy.interpolate import interp1d
from scipy.special import digamma

from structure.units import kB_eV

logger = logging.getLogger(__name__)

try:
    quad_dtype = np.float128  # type: ignore[attr-defined]
except AttributeError:  # pragma: no cover - platforms without float128
    quad_dtype = np.longdouble


class MissingDependencyError(RuntimeError):
    """Raised when optional dependencies required for refinement are absent."""


def _require_mpmath() -> None:
    if mp is None:
        raise MissingDependencyError(
            "mpmath is required for iterative refinement but is not available."
        )


# ---------------------------------------------------------------------------
# Complex trigamma
# ---------------------------------------------------------------------------
#
# psi'(z) is not available for complex arguments in scipy, so this module used
# to call mpmath.polygamma element by element from a Python loop.  With O(1e5)
# probe energies per iteration that loop dominated the entire refinement cost.
#
# The replacement is the textbook recurrence + asymptotic pair, evaluated with
# whole-array numpy operations:
#
#   recurrence   psi'(z) = 1/z^2 + psi'(z+1)          (shift Re z upwards)
#   asymptotic   psi'(z) ~ 1/z + 1/(2 z^2)
#                          + sum_k B_2k / z^(2k+1)    (large |z|)
#
# In this module z = 1/2 + beta/(2 pi) (gamma + i eps) always has Re z >= 1/2,
# so |arg z| < pi/2 and the asymptotic series is well behaved; shifting only
# ever increases Re z and therefore never crosses a pole.
#
# Bernoulli numbers B_2, B_4, ... B_14 as coefficients of z^-(2k+1).
_TRIGAMMA_BERNOULLI = (
    1.0 / 6.0,
    -1.0 / 30.0,
    1.0 / 42.0,
    -1.0 / 30.0,
    5.0 / 66.0,
    -691.0 / 2730.0,
    7.0 / 6.0,
)

# |z| above which the truncated asymptotic series is used.  With terms through
# B_14 the truncation error is ~ |B_16| / |z|^17, i.e. below 1e-16 at |z| = 12.
_TRIGAMMA_ASYMPTOTIC_CUTOFF = 12.0

# Hard cap on recurrence steps; only reached for arguments very close to the
# origin, which cannot occur for the z used here (Re z >= 1/2).
_TRIGAMMA_MAX_SHIFTS = 64


def _trigamma_series(z: np.ndarray) -> np.ndarray:
    """Vectorised complex psi'(z) via recurrence + asymptotic expansion.

    Accurate to ~1e-14 relative for Re(z) > 0, which is far inside the
    tolerance of the sum-rule metric (error_tol >= 1e-4 in practice).
    """

    z = np.asanyarray(z, dtype=np.complex128)
    if np.any(z.real <= 0.0):
        raise ValueError(
            "_trigamma_series requires Re(z) > 0; got Re(z)_min = "
            f"{float(np.min(z.real))}."
        )

    shifted    = z.copy()
    correction = np.zeros_like(shifted)

    # psi'(z) = sum_j 1/(z+j)^2 + psi'(z+n): accumulate until |z+n| is large
    # enough for the asymptotic series.  The number of steps depends only on
    # |z|, so this loop runs at most a handful of times for realistic input.
    for _ in range(_TRIGAMMA_MAX_SHIFTS):
        small = np.abs(shifted) < _TRIGAMMA_ASYMPTOTIC_CUTOFF
        if not small.any():
            break
        correction[small] += 1.0 / shifted[small] ** 2
        shifted[small]    += 1.0
    else:  # pragma: no cover - unreachable for Re(z) >= 1/2
        logger.warning(
            "_trigamma_series hit the recurrence cap (%d shifts); accuracy may "
            "be reduced.", _TRIGAMMA_MAX_SHIFTS,
        )

    inv2 = 1.0 / shifted**2
    # Horner in 1/z^2 for the Bernoulli tail: (B_2 + B_4/z^2 + ...) / z^3
    tail = np.zeros_like(shifted)
    for coeff in reversed(_TRIGAMMA_BERNOULLI):
        tail = tail * inv2 + coeff
    tail *= inv2 / shifted

    return correction + 1.0 / shifted + 0.5 * inv2 + tail


def _trigamma_mpmath(z: np.ndarray) -> np.ndarray:
    """Reference implementation: element-wise mpmath.polygamma(1, z).

    Kept for validation and for the rare case where full mpmath precision is
    wanted.  Orders of magnitude slower than :func:`_trigamma_series`.
    """

    _require_mpmath()
    z = np.asanyarray(z, dtype=np.complex128)
    out = np.empty_like(z, dtype=np.complex128)
    flat_in = z.ravel()
    flat_out = out.ravel()
    for idx, val in enumerate(flat_in):
        flat_out[idx] = complex(mp.polygamma(1, val))
    return out


def trigamma_vec(z: np.ndarray, method: str = "series") -> np.ndarray:
    """Vectorised complex trigamma psi'(z).

    Parameters
    ----------
    z : array_like of complex
        Arguments with Re(z) > 0.
    method : {'series', 'mpmath'}
        ``'series'`` (default) uses the pure-numpy recurrence + asymptotic
        expansion.  ``'mpmath'`` falls back to the element-wise reference
        implementation and requires the optional mpmath dependency.
    """

    if method == "series":
        return _trigamma_series(z)
    if method == "mpmath":
        return _trigamma_mpmath(z)
    raise ValueError(f"unknown trigamma method '{method}'")


def _beta(temperature: float) -> float:
    """Inverse temperature in eV^-1 from a temperature in Kelvin.

    All internal energies (band energies, gamma, mu) are in eV, so beta must
    carry the Boltzmann constant: beta = 1 / (kB * T).  Using 1/T[K] directly
    understates beta by a factor 1/kB ~ 1.16e4, i.e. it corresponds to a
    temperature 1.16e4 times higher than requested.
    """

    if temperature <= 0.0:
        raise ValueError(
            "temperature must be strictly positive (beta = 1/(kB*T) diverges "
            "at T = 0); pass a small finite value instead."
        )
    return 1.0 / (kB_eV * temperature)


def digamma_im(epsilon: np.ndarray, temperature: float, gamma: float) -> np.ndarray:
    """Imaginary part of the digamma-based auxiliary function.

    *temperature* is in Kelvin, *gamma* and *epsilon* in eV.
    """

    beta = _beta(temperature)
    z = 0.5 + (beta / (2.0 * np.pi)) * (gamma + 1j * epsilon)
    return 0.5 - (1.0 / np.pi) * np.imag(digamma(z))


def df_da(epsilon: np.ndarray, temperature: float, gamma: float) -> np.ndarray:
    """d/d(alpha) of the auxiliary function appearing in the error metric.

    *temperature* is in Kelvin, *gamma* and *epsilon* in eV.
    """

    beta = _beta(temperature)
    z = 0.5 + (beta / (2.0 * np.pi)) * (gamma + 1j * epsilon)
    trig = trigamma_vec(z)
    return beta / (2.0 * np.pi**2) * np.real(trig)


@dataclass
class BandData:
    energies:     np.ndarray   # (nkp, nbands)
    k_points:     np.ndarray   # (nkp, 3)  fractional coordinates
    weights:      np.ndarray   # (nkp,)
    multiplicity: np.ndarray   # (nkp,)
    # Per-point Voronoi cell widths along each axis, shape (nkp, 3).
    #
    # For a coarse (ir)reducible mesh produced by ltb run, the initial HDF5
    # does not contain this dataset.  load_band_data then derives the values
    # from the stored nkx/nky/nkz fields: cell_deltas[:,i] = 1/nk_i.  This
    # is correct for both reducible and irreducible meshes because each
    # k-point represents exactly the same Voronoi cell volume (1/nk_i per
    # axis) regardless of symmetry folding.
    #
    # From iteration 1 onward the dataset is written explicitly by
    # write_custom_mesh, so the fallback is never reached.  Each child
    # inherits cell_deltas = parent_cell_deltas / refinement_factor, ensuring
    # the correct cell width is available at every subsequent iteration
    # regardless of which neighbouring points are coarse or fine.
    #
    # NOTE: the nkx/nky/nkz fallback must NOT be used for custom-mesh HDF5
    # files produced by the generator (ltb run with setCustomKmesh), because
    # those always write nkx=nky=nkz=1 — which would give cell_deltas=1.0,
    # i.e. a full BZ width, which is wrong.  write_custom_mesh prevents this
    # by always writing the correct cell_deltas explicitly, and load_band_data
    # raises a hard error if it encounters that combination anyway.
    cell_deltas:  np.ndarray   # (nkp, 3)
    # Where cell_deltas came from: 'file' if read from /.kmesh/cell_deltas
    # (i.e. the input is an already-refined mesh -> continuation), or
    # 'coarse-grid' if derived from nkx/nky/nkz (fresh coarse input).
    cell_deltas_source: str = "file"


def load_band_data(h5file: "h5py.File") -> BandData:
    """Load band energies and k-mesh data from an HDF5 handle.

    cell_deltas is read from the ``/.kmesh/cell_deltas`` dataset when present.
    This dataset is written by :func:`write_custom_mesh` on every iteration
    and patched into every generator output by the refinement loop; its
    presence marks the file as an already-(partially-)refined mesh, from
    which refinement can be *continued* (``cell_deltas_source == 'file'``).

    For an initial coarse HDF5 (produced by ``ltb run``, where cell_deltas is
    absent), the values are derived from the stored ``nkx``/``nky``/``nkz``
    fields as ``1/nk_i`` (``cell_deltas_source == 'coarse-grid'``).

    A custom-mesh file (nkx=nky=nkz=1) without cell_deltas is rejected with a
    ValueError: its per-point cell widths are unknown and the fallback would
    silently produce zero-width cells, i.e. a no-op refinement.
    """

    energies = h5file["/energies"][:].astype(quad_dtype)
    if energies.ndim != 2 or energies.shape[1] < 1:
        raise ValueError("Energy dataset must be 2D with at least one band.")

    k_group      = h5file["/.kmesh"]
    k_points     = k_group["points"][:].astype(quad_dtype)
    weights      = k_group["weights"][:].astype(quad_dtype)
    multiplicity = k_group["multiplicity"][:].astype(quad_dtype)

    if "cell_deltas" in k_group:
        # Present in every generator output from iteration 1 onward (patched
        # by the refinement loop).  Its presence therefore identifies the
        # input as an already-(partially-)refined mesh -> continuation.
        cell_deltas = k_group["cell_deltas"][:].astype(quad_dtype)
        if cell_deltas.shape != (k_points.shape[0], 3):
            raise ValueError(
                "/.kmesh/cell_deltas has shape {} but the mesh has {} k-points; "
                "the file appears corrupted or was assembled from mismatched "
                "refinement stages.".format(cell_deltas.shape, k_points.shape[0])
            )
        source = "file"
    else:
        nkx = int(k_group["nkx"][()])
        nky = int(k_group["nky"][()])
        nkz = int(k_group["nkz"][()])
        if nkx * nky * nkz == 1 and k_points.shape[0] > 1:
            # Custom-mesh convention: generator outputs always write
            # nkx=nky=nkz=1.  Deriving cell widths from these would give
            # zero-width cells and the refinement would silently do nothing.
            # This combination means the file is a refined/custom mesh whose
            # cell_deltas dataset is missing (produced by a pipeline version
            # predating the cell_deltas patch step, or stripped by an
            # external tool).  Refuse loudly instead of no-op'ing.
            raise ValueError(
                "Input HDF5 is a custom-mesh file (nkx=nky=nkz=1, {} k-points) "
                "but contains no /.kmesh/cell_deltas dataset. Refinement cannot "
                "be (re)started from it: per-point cell widths are unknown. "
                "Either restart from the original coarse mesh, or restore the "
                "cell_deltas dataset from the corresponding custom_mesh_iter_N "
                "file of the run that produced it.".format(k_points.shape[0])
            )
        # Initial coarse HDF5 from ltb run: nkx/nky/nkz are the true coarse
        # grid dimensions.  Derive uniform cell widths from them.
        cell_deltas = np.tile(
            np.array(
                [
                    quad_dtype(1.0) / quad_dtype(nkx) if nkx > 1 else quad_dtype(0.0),
                    quad_dtype(1.0) / quad_dtype(nky) if nky > 1 else quad_dtype(0.0),
                    quad_dtype(1.0) / quad_dtype(nkz) if nkz > 1 else quad_dtype(0.0),
                ],
                dtype=quad_dtype,
            ),
            (k_points.shape[0], 1),
        )
        source = "coarse-grid"

    return BandData(
        energies=energies,
        k_points=k_points,
        weights=weights,
        multiplicity=multiplicity,
        cell_deltas=cell_deltas,
        cell_deltas_source=source,
    )


def build_band_axis(bands: np.ndarray, max_bands: "int | None" = None) -> np.ndarray:
    """Construct a monotonic energy axis from the sampled band energies.

    Parameters
    ----------
    bands : ndarray (nkp, nbands)
        Band energies of the current mesh.
    max_bands : int or None
        If None (default) *all* bands enter the axis.  This is the correct
        choice for the df/da sum-rule metric: the bands crossing mu are not
        in general the lowest ones, and restricting the axis to a subset both
        removes the energies that actually resolve the peak at mu and
        truncates the integration range, which adds a spurious error floor.
        An explicit integer restricts the axis to the ``max_bands`` lowest
        bands (legacy behaviour, kept for callers that need it).
    """

    if max_bands is None:
        subset = bands
    else:
        take = min(max_bands, bands.shape[1])
        subset = bands[:, :take]
    combined = np.sort(np.unique(subset.reshape(-1)))
    return combined.astype(quad_dtype)


def compute_error(band: np.ndarray, temperature: float, gamma: float) -> float:
    """Return |1 - int df_da d epsilon| for the supplied band grid."""

    dfda = df_da(band, temperature, gamma)
    # np.trapz was removed in numpy 2.x in favour of np.trapezoid
    _trapz = getattr(np, "trapezoid", None) or np.trapz
    integral = _trapz(dfda, band)
    return float(np.abs(1.0 - integral))


def kernel_width(temperature: float, gamma: float) -> float:
    """Half-width [eV] of the df/da peak at mu: max(kB*T, gamma).

    This is the only energy scale the sum-rule metric can resolve.  A k-mesh
    that does not place several distinct band energies inside this width
    cannot satisfy the sum rule regardless of how hotspots are selected, and
    the hotspot window must never shrink below it.
    """

    return max(kB_eV * float(temperature), float(gamma))


def _build_probe_grid(
    band_min:           float,
    band_max:           float,
    chemical_potential: float,
    energy_window:      float,
    width:              float,
    log_points:         int,
    uniform_points:     int,
) -> np.ndarray:
    """Probe energies: uniform over the bandwidth plus log-clustered at mu.

    The clustering is built from *offsets* relative to mu that are
    geometrically spaced between a small fraction of the kernel width and
    ``energy_window``.  The previous version spaced the offsets as
    ``geomspace(|mu|, |mu| + window)``, which is only genuinely logarithmic
    for mu ~ 0: for any other chemical potential the ratio
    (|mu| + window)/|mu| is O(1) and the "log" grid degenerates into a second
    uniform grid, leaving the peak core unresolved.
    """

    uniform_band = np.linspace(band_min, band_max, num=uniform_points,
                               endpoint=True, dtype=float)

    inner = max(1e-4 * width, 1e-12)
    outer = max(float(energy_window), 10.0 * inner)
    offsets = np.geomspace(inner, outer, num=log_points)

    return np.sort(np.unique(np.concatenate((
        uniform_band,
        chemical_potential + offsets,
        chemical_potential - offsets,
        np.array([chemical_potential], dtype=float),
    ))))


# Cache for the probe grid and the analytic kernel evaluated on it.  Both
# depend only on (band range, mu, window, T, gamma, point counts) -- none of
# which change during a refinement run except the band range, which is
# usually constant as well.  Recomputing them every iteration was pure waste.
_PROBE_CACHE: "dict[tuple, tuple[np.ndarray, np.ndarray]]" = {}
_PROBE_CACHE_MAXSIZE = 8


def _probe_and_analytic(
    band_min:           float,
    band_max:           float,
    temperature:        float,
    gamma:              float,
    chemical_potential: float,
    energy_window:      float,
    log_points:         int,
    uniform_points:     int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (probe_band, analytic df/da on it), memoised."""

    key = (band_min, band_max, temperature, gamma, chemical_potential,
           energy_window, log_points, uniform_points)
    hit = _PROBE_CACHE.get(key)
    if hit is not None:
        return hit

    width      = kernel_width(temperature, gamma)
    probe_band = _build_probe_grid(band_min, band_max, chemical_potential,
                                   energy_window, width,
                                   log_points, uniform_points)
    analytic   = df_da(probe_band, temperature, gamma)

    if len(_PROBE_CACHE) >= _PROBE_CACHE_MAXSIZE:
        _PROBE_CACHE.pop(next(iter(_PROBE_CACHE)))
    _PROBE_CACHE[key] = (probe_band, analytic)
    return probe_band, analytic


def clear_probe_cache() -> None:
    """Drop the memoised probe grids (useful in long-running processes)."""

    _PROBE_CACHE.clear()


def detect_refinement_scale(
    band: np.ndarray,
    temperature: float,
    gamma: float,
    chemical_potential: float,
    energy_window: float = 0.1,
    log_points: int = 12_000,
    uniform_points: int = 4_000,
    defect_rtol: float = 1e-3,
    width_factor: float = 2.0,
) -> Tuple[np.ndarray, float]:
    """Identify undersampled energies and derive an adaptive hotspot window.

    The current band axis is used as the node set of a linear interpolant of
    df/da.  Where that interpolant deviates from the analytic kernel, the mesh
    is undersampled in energy.  The *reach* of that deviation around mu sets
    the energy window within which k-points are flagged for subdivision.

    Window policy ("narrow only")
    -----------------------------
    ``tolerance = clip(reach, floor, energy_window)`` with
    ``reach = max|significant - mu|`` taken over the significant energies that
    lie *inside* ``energy_window`` (band-edge van Hove defects are always
    present and would otherwise pin the window at its ceiling), and
    ``floor = min(energy_window, width_factor * kernel_width)``.

    The user-supplied ``energy_window`` is a hard ceiling, so the window can
    only ever shrink relative to it as the mesh improves.  This keeps the
    adaptivity that the defect analysis provides while making the historical
    failure mode structurally impossible: an earlier version scaled the
    window as ``3 * max|significant|``, which for a Dirac point at mu spanned
    the whole bandwidth and flagged every k-point.  The floor prevents the
    opposite failure -- a window narrower than the kernel itself, which would
    stop refining points that still carry weight in the integral.

    Parameters
    ----------
    band : ndarray
        Sorted unique band energies of the current mesh.
    temperature, gamma : float
        Kelvin and eV; the worst-case corner of the planned sweep.
    chemical_potential : float
        mu [eV].
    energy_window : float
        Ceiling for the hotspot window [eV].
    log_points, uniform_points : int
        Size of the probe grid.  Reduced from (110000, 10000): the metric only
        needs max|defect| and its half-level set, which are resolved to well
        under a percent at these sizes, and the analytic kernel on the grid is
        cached across iterations anyway.
    defect_rtol : float
        Relative defect below which the band axis is considered to resolve the
        kernel everywhere.  Returning ``tolerance = 0.0`` then signals the
        caller to stop: no amount of extra k-points near mu can reduce the
        sum-rule error further, which is the plateau situation.
    width_factor : float
        Multiplier of the kernel width used for the window floor.

    Returns
    -------
    significant : ndarray
        Probe energies whose interpolation defect is at least half the maximum.
    tolerance : float
        Hotspot half-window [eV] for :func:`refine_kmesh`, or 0.0 if the band
        axis already resolves the kernel to within ``defect_rtol``.
    """

    dfda = df_da(band, temperature, gamma)

    band_min = float(np.min(band))
    band_max = float(np.max(band))

    probe_band, analytic = _probe_and_analytic(
        band_min, band_max, float(temperature), float(gamma),
        float(chemical_potential), float(energy_window),
        log_points, uniform_points,
    )

    interpolator = interp1d(band, dfda, kind="linear",
                            bounds_error=False, fill_value="extrapolate")
    interpolated = interpolator(probe_band)
    delta_arr    = np.abs(analytic - interpolated)

    if delta_arr.size == 0:
        return probe_band, 0.0

    max_defect = float(np.max(delta_arr))
    max_kernel = float(np.max(np.abs(analytic)))

    # Converged in energy space: the sampled energies already reproduce the
    # kernel everywhere, so further k-points near mu cannot help.
    if max_kernel > 0.0 and max_defect <= defect_rtol * max_kernel:
        logger.info(
            "Band axis resolves df/da to %.3e relative (<= defect_rtol=%.3e); "
            "no undersampled energies remain.",
            max_defect / max_kernel, defect_rtol,
        )
        return probe_band, 0.0

    significant = probe_band[delta_arr >= 0.5 * max_defect]

    if significant.size == 0:
        return probe_band, 0.0

    width = kernel_width(temperature, gamma)
    floor = min(float(energy_window), width_factor * width)

    # Only defects *inside* the user window may set the window: band-edge
    # (van Hove) defects are always present and far from mu, and letting them
    # enter the maximum would pin the window open at the ceiling forever.
    offsets   = np.abs(significant - chemical_potential)
    in_window = offsets <= float(energy_window)
    reach     = float(np.max(offsets[in_window])) if in_window.any() else 0.0

    tolerance = min(float(energy_window), max(floor, reach))
    tolerance = max(tolerance, 1e-6)

    logger.debug(
        "detect_refinement_scale: max defect %.3e (%.2f%% of peak), "
        "reach %.4e eV, floor %.4e eV, window %.4e eV -> tolerance %.4e eV.",
        max_defect, 100.0 * max_defect / max_kernel if max_kernel else float('nan'),
        reach, floor, energy_window, tolerance,
    )

    return significant.astype(quad_dtype), tolerance


def _merge_duplicates_with_tolerance(
    points:      np.ndarray,
    weights:     np.ndarray,
    cell_deltas: np.ndarray,
    tol:         float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if points.size == 0:
        return points, weights, cell_deltas

    scaled     = np.round(points / tol) * tol
    dtype      = scaled.dtype
    structured = scaled.view([(f"f{i}", dtype) for i in range(scaled.shape[1])])
    _, idx, inv = np.unique(structured, return_index=True, return_inverse=True)
    # numpy >= 2.0 returns inverse indices with the shape of the input array
    # (here (M, 1) from the structured view) instead of flattened -> ravel
    inv = np.asarray(inv).ravel()
    merged_points  = points[idx]
    merged_weights = np.zeros(len(idx), dtype=weights.dtype)
    np.add.at(merged_weights, inv, weights)
    # For cell_deltas keep the first occurrence: all duplicates are the same
    # physical point so they carry identical deltas by construction.
    merged_deltas = cell_deltas[idx]
    return merged_points, merged_weights, merged_deltas


def refine_kmesh(
    band_data:         BandData,
    target_energy:     float,
    tolerance:         float,
    refinement_factor: int   = 3,
    uniqueness_tol:    float = 1e-12,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return refined (points, weights, cell_deltas) preserving total weight.

    Each hotspot k-point is replaced by ``refinement_factor^3`` children
    placed symmetrically around the parent:

        k_child = k_parent + (i - (n-1)/2) * cell_deltas[parent] / n

    for i = 0 … n-1 along each active axis.  For odd n the central child
    coincides with the parent, so the parent is always included in the
    refined mesh.  Even refinement_factor values are incremented by 1.

    Each child inherits ``cell_deltas = parent_cell_deltas / n``, which is
    written to HDF5 by :func:`write_custom_mesh` and read back by
    :func:`load_band_data` on the next iteration.  This guarantees that the
    correct Voronoi cell width is used at every iteration regardless of the
    local mesh density around the point being refined.

    Weights are split *per parent*: each child of a parent carries
    ``w_parent / n_children(parent)``.  A parent's cell is exactly tiled by
    its own children, so this is the only split that conserves weight
    locally.  (Redistributing the pooled removed weight uniformly over all
    children of all parents conserves the total as well, but transfers weight
    between parents whenever their weights differ -- which is the generic
    case for an irreducible input mesh, where the weight encodes the star
    multiplicity.)

    Returns
    -------
    merged_points  : ndarray (M, 3)
    merged_weights : ndarray (M,)
    merged_deltas  : ndarray (M, 3)  — must be persisted via write_custom_mesh
    """

    if refinement_factor % 2 == 0:
        refinement_factor += 1
        logger.warning(
            "refinement_factor must be odd for symmetric cell-centred "
            "subdivision.  Incrementing to %d.", refinement_factor
        )

    energies    = band_data.energies
    k_points    = band_data.k_points
    weights     = band_data.weights
    cell_deltas = band_data.cell_deltas

    tol    = quad_dtype(tolerance)
    target = quad_dtype(target_energy)

    # Explicit absolute window.  np.isclose would add its default relative
    # term rtol*|target| = 1e-5*|mu|, silently widening the window by an
    # amount that depends on the (arbitrary) energy zero of the model.
    mask    = np.any(np.abs(energies - target) <= tol, axis=1)
    indices = np.where(mask)[0]

    if indices.size == 0:
        logger.warning(
            "No k-points within tolerance %.3e of target energy %.3e.",
            tolerance, target_energy,
        )
        return k_points.copy(), weights.copy(), cell_deltas.copy()

    keep_mask                  = np.ones(k_points.shape[0], dtype=bool)
    removed_weights            = []
    refined_points_collection  = []
    refined_deltas_collection  = []
    refined_weights_collection = []

    one     = quad_dtype(1.0)
    n       = refinement_factor
    offsets = np.arange(n, dtype=quad_dtype) - quad_dtype(n - 1) / 2

    for idx in indices:
        keep_mask[idx] = False
        removed_weights.append(weights[idx])

        k_pt        = k_points[idx]
        deltas      = cell_deltas[idx]
        child_deltas = deltas / quad_dtype(n)

        refined_axes = [
            k_pt[axis] + offsets * child_deltas[axis]
            if deltas[axis] > 0 else np.array([k_pt[axis]], dtype=quad_dtype)
            for axis in range(3)
        ]

        mesh    = np.meshgrid(*refined_axes, indexing="ij")
        refined = np.column_stack([ax.ravel() for ax in mesh]).astype(quad_dtype)
        refined %= one
        refined[refined >= (one - quad_dtype(1e-15))] = quad_dtype(0.0)
        refined_points_collection.append(refined)

        n_children = refined.shape[0]
        refined_deltas_collection.append(np.tile(child_deltas, (n_children, 1)))
        # per-parent weight split: the parent cell is exactly tiled by these
        # n_children sub-cells, so each carries 1/n_children of its weight.
        refined_weights_collection.append(
            np.full(n_children, weights[idx] / quad_dtype(n_children),
                    dtype=quad_dtype)
        )

    k_points_filtered = k_points[keep_mask]
    weights_filtered  = weights[keep_mask]
    deltas_filtered   = cell_deltas[keep_mask]

    total_removed = np.sum(removed_weights) if removed_weights else quad_dtype(0.0)

    if refined_points_collection:
        refined_points  = np.concatenate(refined_points_collection,  axis=0)
        refined_deltas  = np.concatenate(refined_deltas_collection,  axis=0)
        refined_weights = np.concatenate(refined_weights_collection, axis=0)
        if not np.allclose(np.sum(refined_weights).astype(float),
                           np.float64(total_removed), rtol=1e-12, atol=1e-14):
            logger.warning(
                "Per-parent weight split does not reproduce the removed "
                "weight: removed=%s redistributed=%s",
                float(total_removed), float(np.sum(refined_weights)),
            )
    else:
        refined_points  = np.empty((0, 3), dtype=quad_dtype)
        refined_deltas  = np.empty((0, 3), dtype=quad_dtype)
        refined_weights = np.empty((0,),   dtype=quad_dtype)

    merged_points  = np.concatenate((k_points_filtered, refined_points),  axis=0)
    merged_weights = np.concatenate((weights_filtered,  refined_weights), axis=0)
    merged_deltas  = np.concatenate((deltas_filtered,   refined_deltas),  axis=0)

    merged_points, merged_weights, merged_deltas = _merge_duplicates_with_tolerance(
        merged_points, merged_weights, merged_deltas, tol=quad_dtype(uniqueness_tol)
    )

    return merged_points, merged_weights, merged_deltas


def write_custom_mesh(
    path:        str,
    points:      np.ndarray,
    weights:     np.ndarray,
    cell_deltas: np.ndarray,
) -> None:
    """Write a custom k-mesh HDF5 file including per-point cell widths.

    ``cell_deltas`` must always be written so that :func:`load_band_data`
    never falls back to the nkx/nky/nkz path on subsequent iterations.
    Generator-produced HDF5 files always contain nkx=nky=nkz=1 (custom mesh
    convention), so the fallback would give cell_deltas=1.0 — the full BZ
    width — which is incorrect.
    """
    import h5py

    with h5py.File(path, "w") as h5file:
        km_group = h5file.create_group(".kmesh")
        km_group.create_dataset("points",      data=points.astype(float))
        km_group.create_dataset("weights",     data=weights.astype(float))
        km_group.create_dataset("cell_deltas", data=cell_deltas.astype(float))
