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


# Gauss-Legendre nodes for the per-panel reference integral of df/da.
# Five nodes is converged: on a 714k-panel graphene mesh the L1 metric is
# 0.006653 with 5 nodes against 0.006652 with 20, i.e. four significant
# figures, at a fifth of the cost (2.6 s vs 14 s).
_GL_NODES = 5
_GL_X, _GL_W = np.polynomial.legendre.leggauss(_GL_NODES)

# Panels processed per chunk.  Each chunk holds an (n, _GL_NODES) float64
# temporary, i.e. ~80 MB at 2e6 panels.
_PANEL_CHUNK = 2_000_000


def panel_defects(
    band:        np.ndarray,
    temperature: float,
    gamma:       float,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Per-panel trapezoid defect of df/da on the band axis.

    For each panel [e_i, e_i+1] of the sorted band axis this returns the
    SIGNED defect

        d_i = trapezoid_i - exact_i,   exact_i = int_{e_i}^{e_i+1} df/da de

    with ``exact_i`` evaluated by ``_GL_NODES``-point Gauss-Legendre.  The
    analytic kernel is available in closed form, so there is no reason to
    estimate the defect by an embedded (bisection) rule: that estimator
    collapses exactly where it is needed most, on a coarse mesh whose panels
    are far wider than the kernel, because neither the trapezoid nor its
    bisection resolves the peak.  Measured on the graphene 48x48 start it
    reported 20.3 against a true defect of 39.7.

    Returns
    -------
    midpoints : ndarray (n_panels,)
        Panel centres, used to locate the defect in energy.
    defects : ndarray (n_panels,)
        Signed defects d_i.
    exact_total : float
        sum_i exact_i, i.e. the exact integral of df/da over the SAMPLED
        energy range.  1 - exact_total is the truncation error of that range.
    """

    band = np.asarray(band, dtype=np.float64).ravel()
    if band.size < 2:
        return (np.empty(0), np.empty(0), 0.0)

    lo, hi = band[:-1], band[1:]
    half   = 0.5 * (hi - lo)
    mid    = 0.5 * (lo + hi)

    node_f = df_da(band, temperature, gamma)
    trap   = half * (node_f[:-1] + node_f[1:])

    exact = np.empty_like(trap)
    for start in range(0, trap.size, _PANEL_CHUNK):
        stop  = min(start + _PANEL_CHUNK, trap.size)
        m, hw = mid[start:stop, None], half[start:stop, None]
        nodes = m + hw * _GL_X[None, :]
        exact[start:stop] = half[start:stop] * np.sum(
            df_da(nodes, temperature, gamma) * _GL_W[None, :], axis=1
        )

    return mid, trap - exact, float(np.sum(exact))


def compute_error(band: np.ndarray, temperature: float, gamma: float) -> float:
    """Quadrature defect of the df/da sum rule on the supplied band axis.

        error = sum_i |d_i|  +  |1 - sum_i exact_i|
                (discretisation)    (truncation of the sampled range)

    The first term replaces the historical ``|1 - int df/da|``, which was the
    absolute value of the SIGNED sum of the same defects.  df/da is a peak at
    mu, concave over its top and convex in both tails, so the trapezoid rule
    underestimates near mu and overestimates in the tails: the two
    contributions ALWAYS carry opposite signs, in every system.  Whenever they
    happen to be of comparable magnitude the signed sum can dip on one iterate
    through accidental cancellation and rise again on the next, strictly finer
    one -- refinement only ever adds k-points, so a rising signed metric never
    means a worse mesh.  Graphene 48x48x1 at gamma = 1 meV showed exactly this:

        nE      signed    this metric
        475    0.04480      0.04501
        587    0.02301      0.04242     <- signed dips (cancellation)
        1435   0.02966      0.03135     <- signed rebounds, L1 keeps falling
        714051 0.00644      0.00665

    Because |sum d_i| <= sum |d_i|, the new value is never smaller than the old
    one: a mesh that passes this test would also have passed the old one, so
    nothing previously converged becomes unconverged.  The converse is the
    price -- runs that used to pass through cancellation now need more
    iterations.

    NOTE this is an energy-axis criterion only: it never sees the k-point
    weights.  Two meshes with identical energy sets score identically however
    differently they sample the Brillouin zone.  Treat it as "is the kernel
    resolved in energy", not as "is sigma converged".
    """

    _, defects, exact_total = panel_defects(band, temperature, gamma)
    if defects.size == 0:
        return float("inf")
    return float(np.sum(np.abs(defects)) + abs(1.0 - exact_total))


def kernel_width(temperature: float, gamma: float) -> float:
    """Half-width [eV] of the df/da peak at mu: max(kB*T, gamma).

    This is the only energy scale the sum-rule metric can resolve.  A k-mesh
    that does not place several distinct band energies inside this width
    cannot satisfy the sum rule regardless of how hotspots are selected, and
    the hotspot window must never shrink below it.
    """

    return max(kB_eV * float(temperature), float(gamma))


# Default ratio between successive kernel widths of the refinement cascade.
# Three matches the default refinement_factor: one ladder step then costs
# roughly one subdivision of the shell that the step exposes.
_LADDER_RATIO_DEFAULT = 3.0
_LADDER_MAX_STAGES_DEFAULT = 6


def kernel_ladder(
    T_min:      float,
    gamma_min:  float,
    T_max:      "float | None" = None,
    gamma_max:  "float | None" = None,
    ratio:      float = _LADDER_RATIO_DEFAULT,
    max_stages: int   = _LADDER_MAX_STAGES_DEFAULT,
) -> "list[tuple[float, float]]":
    """Kernel-width ladder for the refinement cascade, WIDEST FIRST.

    Returns ``[(T_0, gamma_0), ..., (T_min, gamma_min)]`` with strictly
    decreasing ``kernel_width``.  The last entry is always exactly
    ``(T_min, gamma_min)``, so a single-stage ladder (the default, when no
    maxima are given) reproduces the pre-cascade behaviour exactly.

    Why widest first
    ----------------
    The wide stage places k-points across the whole ``+-W_max`` shell while
    the shell is still cheap to cover, and those points become the parents
    that the narrow stages subdivide.  The alternative -- taking ``max()`` of
    the metric over the ladder inside a single loop -- does NOT do this: on a
    coarse mesh the narrow-kernel defect dwarfs the wide one (39.7 against
    much less on the graphene 48x48 start), so the maximum is always the
    narrow entry and the loop degenerates into today's narrow-first
    behaviour.  The ordering has to be explicit.

    Spacing
    -------
    Stages are spaced geometrically in ``W = max(kB*T, gamma)``, i.e.
    uniformly in ``log W`` -- the two descriptions are the same thing.  The
    *ratio* is pinned rather than the stage count, so a wider span costs more
    stages instead of larger jumps per stage; ``max_stages`` caps the total.
    ``T`` and ``gamma`` are interpolated independently on the same normalised
    parameter, so a ladder that moves only one of them leaves the other fixed.

    Parameters
    ----------
    T_min, gamma_min : float
        The narrow (target) corner: Kelvin and eV.
    T_max, gamma_max : float or None
        The wide corner.  ``None`` means "same as the minimum", which yields a
        one-entry ladder.  Values below the corresponding minimum are raised
        to it with a warning: the ladder must end at the target corner.
    ratio : float
        Ratio between the widths of successive stages (> 1).
    max_stages : int
        Hard cap on the number of stages.

    Returns
    -------
    list of (temperature, gamma) tuples, widest kernel first.
    """

    T_min     = float(T_min)
    gamma_min = float(gamma_min)
    T_max     = T_min     if T_max     is None else float(T_max)
    gamma_max = gamma_min if gamma_max is None else float(gamma_max)

    if T_max < T_min:
        logger.warning(
            "T_max = %.6g K is below T_min = %.6g K; the cascade must end at "
            "the target corner, so T_max is raised to T_min.", T_max, T_min,
        )
        T_max = T_min
    if gamma_max < gamma_min:
        logger.warning(
            "gamma_max = %.4e eV is below gamma_min = %.4e eV; the cascade "
            "must end at the target corner, so gamma_max is raised to "
            "gamma_min.", gamma_max, gamma_min,
        )
        gamma_max = gamma_min

    if ratio <= 1.0:
        raise ValueError("ladder ratio must be greater than 1, got %r" % ratio)
    if max_stages < 1:
        raise ValueError("max_stages must be at least 1, got %r" % max_stages)

    W_min = kernel_width(T_min, gamma_min)
    W_max = kernel_width(T_max, gamma_max)

    # Nothing to cascade over: the two corners resolve the same energy scale.
    if W_min <= 0.0 or W_max <= W_min * (1.0 + 1e-12):
        return [(T_min, gamma_min)]

    n_stages = int(np.ceil(np.log(W_max / W_min) / np.log(ratio))) + 1
    if n_stages > max_stages:
        logger.info(
            "Kernel ladder from %.4e to %.4e eV at ratio %.3g would need %d "
            "stages; capped at max_stages = %d (effective ratio %.3g).",
            W_max, W_min, ratio, n_stages, max_stages,
            (W_max / W_min) ** (1.0 / max(max_stages - 1, 1)),
        )
        n_stages = max_stages

    # The ladder is geometric in the WIDTH, which is the only energy scale the
    # metric resolves (see kernel_width).  Interpolating T and gamma
    # independently does NOT achieve that: max(kB*T, gamma) switches which of
    # the two controls the width part-way down the ladder, and the widths then
    # come out unevenly spaced (measured 5.57, 2.15, 2.15 for a nominal ratio
    # of 3).  So the widths are laid out first, and each is realised on
    # whichever channel is free to carry it:
    #
    #   gamma_j = clip(W_j,         gamma_min, gamma_max)
    #   T_j     = clip(W_j / kB,    T_min,     T_max)
    #
    # each clipped to its own requested range and held at its minimum when
    # that channel was not asked to vary.  Because W_max is by construction
    # attained at the wide corner, at least one channel reaches W_j for every
    # rung, so max(kB*T_j, gamma_j) == W_j exactly; and a gamma-only cascade
    # leaves the temperature alone, as a user who moved only gamma expects.
    widths = np.geomspace(W_max, W_min, n_stages)

    T_varies     = T_max     > T_min
    gamma_varies = gamma_max > gamma_min

    ladder = []
    for w in widths:
        t = (min(max(float(w) / kB_eV, T_min), T_max)) if T_varies else T_min
        g = (min(max(float(w), gamma_min), gamma_max)) if gamma_varies else gamma_min
        ladder.append((float(t), float(g)))

    # Guard the endpoints against float drift in the geomspace/clip round trip.
    ladder[0]  = (T_max, gamma_max)
    ladder[-1] = (T_min, gamma_min)
    return ladder


def stage_tolerances(
    ladder:    "list[tuple[float, float]]",
    error_tol: float,
    exponent:  float = 0.0,
) -> "list[float]":
    """Per-stage error tolerance for a cascade ladder.

        tol_j = error_tol * (W_j / W_min) ** exponent

    ``exponent = 0`` (the default) applies the SAME tolerance at every kernel
    width.  That is not merely the conservative choice, it is the physically
    natural one: the per-panel trapezoid defect of a kernel of width ``W``
    sampled with panel width ``h`` is

        eps ~ (|g''|/12) * (h/W)**2                     (see kernel_width)

    -- a function of ``h/W`` alone.  Holding ``eps`` fixed across the ladder
    therefore holds ``h/W`` fixed, i.e. produces a mesh whose local energy
    spacing is proportional to its distance from mu.  That self-similar
    grading is exactly what a peaked kernel with power-law tails wants.

    It is also affordable.  For a ``(d-1)``-dimensional Fermi surface in a
    ``d``-dimensional zone the shell ``|e - mu| < W`` has volume ``~ W`` and
    is covered by cells of side ``h/v``, so a stage costs
    ``N_j ~ W_j / h_j**d ~ W_j**(1-d)`` at ``exponent = 0``: the wide stages
    are CHEAPER than the narrow one (``1/W`` in 2D, ``1/W**2`` in 3D) and the
    cascade costs barely more than the narrow stage alone.  Only for a point
    node (graphene at the Dirac point, a 0-dimensional Fermi surface) is
    ``N_j`` constant per stage, and the cascade then costs ``n_stages`` times
    the narrow-only run -- on a few-hundred-point mesh, nothing.

    ``exponent > 0`` loosens the wide stages geometrically along with the
    width and is offered for cases where the wide corner is only meant to
    bracket.  It is not free: points not placed in the ``W_max`` shell while
    that stage is running are never recovered, because at ``W_min`` those
    panels carry negligible defect and are never selected by
    :func:`select_refinement_panels`.  The final mesh is then converged at
    ``W_min`` and loose by ``(W_max/W_min)**exponent`` at ``W_max``.
    """

    if not ladder:
        return []
    widths = [kernel_width(t, g) for t, g in ladder]
    w_min  = widths[-1]
    if w_min <= 0.0 or exponent == 0.0:
        return [float(error_tol)] * len(ladder)
    return [float(error_tol) * (w / w_min) ** float(exponent) for w in widths]


def _symmetry_axis_orbits(symop: np.ndarray, dims: np.ndarray) -> "list[list[int]]":
    """Partition the three reciprocal axes into orbits under the point group.

    Axes that a symmetry operation maps into one another MUST carry the same
    number of divisions, or the resulting grid is not invariant under the
    group and the irreducible reduction silently falls back to a smaller
    subgroup.  ``symop`` is the stack of integer rotation matrices stored in
    ``/.unitcell/symop``; an operation connects axes i and j when its matrix
    has a non-zero entry coupling them.

    Inactive dimensions (``dims[i] is False``) are returned as singleton
    orbits and are never merged with active ones.
    """

    parent = list(range(3))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    ops = np.asarray(symop)
    if ops.size:
        ops = ops.reshape(-1, 3, 3)
        for op in ops:
            for i in range(3):
                for j in range(3):
                    if i != j and abs(op[i, j]) > 1e-8 and dims[i] and dims[j]:
                        union(i, j)

    orbits: "dict[int, list[int]]" = {}
    for ax in range(3):
        orbits.setdefault(find(ax), []).append(ax)
    return list(orbits.values())


def axis_anisotropy(
    energies:    np.ndarray,
    moments:     np.ndarray,
    weights:     np.ndarray,
    kvec:        np.ndarray,
    dims:        np.ndarray,
    temperature: float,
    gamma:       float,
    symop:       "np.ndarray | None" = None,
) -> np.ndarray:
    """RMS energy variation per unit fractional step along each reciprocal axis.

    The natural number of divisions along axis ``i`` is proportional to the
    returned ``d[i]``.  Minimising the total discretisation error, which goes
    as ``sum_j (d_j / W)**2``, at fixed cost ``N = prod_j n_j`` gives
    ``n_j**3 ~ c_j**2 * n_j`` and hence ``n_j ~ d_j`` -- LINEAR, with no
    additional velocity weighting, which would over-steepen the ratio.

    The estimate uses the band-diagonal velocity tensor rather than finite
    differences between neighbouring k-points:

        d_i**2 = sum_kn w_kn (b_i . M(k,n) . b_i) / sum_kn w_kn
        w_kn   = weights_k * df/da(E_n(k) - mu ; W)

    which has four properties finite differences do not:

    * No neighbour structure is needed, so it works directly on an
      IRREDUCIBLE wedge, on a refined mesh, or on any custom point set.
    * Band selection is automatic.  The kernel weight at the WIDE corner
      suppresses bands far from mu, so each band contributes in proportion to
      what it actually contributes to transport; no band window is needed,
      and no band-ordering problem arises because nothing is differenced
      across k.
    * Non-orthogonal lattices are handled, because the full tensor is
      contracted with the reciprocal vectors instead of picking out diagonal
      Cartesian components.  This is why graphene returns exactly 1:1 on
      hexagonal axes.
    * It is gauge-safe at degeneracies.  Per KI-11 the band-diagonal moments
      are gauge-arbitrary within a degenerate multiplet and only their sum is
      well defined -- but degenerate partners share an energy, hence an
      identical kernel weight, so the weighted sum is gauge-invariant.

    Note on accuracy: the estimate is a good predictor but should be treated
    as a LOWER bound.  On an orthorhombic model with t_x/t_y = 5 it returns
    6.97, which sits right at the onset of the accuracy plateau: at matched
    cost 1:1 gives +12.3%, 2:1 +2.1%, 4.6:1 +0.50%, 8:1 +0.004%, and the
    plateau extends to at least 32:1.  Over-shooting is therefore nearly free
    while under-shooting is not, so :func:`suggest_mesh_ratio` rounds up.
    The estimate is settled by about nk = 32 per axis; a coarser trial mesh
    under-estimates it (6.09 at nk = 16), again erring in the safe direction.

    Parameters
    ----------
    energies : (nkp, nbands) array, already referenced to mu.
    moments  : (nkp, nbands, ncomp) band-diagonal moments; ncomp is 3
        (xx, yy, zz) for orthogonal cells and 6 (xx, yy, zz, xy, xz, yz)
        otherwise.
    weights  : (nkp,) k-point weights.
    kvec     : (3, 3) reciprocal lattice vectors, one per row.
    dims     : (3,) boolean mask of active dimensions.
    temperature, gamma : the WIDE corner of the planned sweep.
    symop    : optional (nsym, 3, 3) rotations; when given, axes in the same
        symmetry orbit are averaged so that symmetry-equivalent directions
        come out exactly equal rather than merely nearly so.

    Returns
    -------
    (3,) array; entries for inactive dimensions are 0.
    """

    energies = np.asarray(energies, dtype=float)
    moments  = np.real(np.asarray(moments))
    weights  = np.asarray(weights, dtype=float)
    kvec     = np.asarray(kvec, dtype=float).reshape(3, 3)
    dims     = np.asarray(dims).astype(bool).ravel()[:3]

    omega = weights[:, None] * df_da(energies, temperature, gamma)
    total = float(np.sum(omega))
    if not np.isfinite(total) or total <= 0.0:
        logger.warning(
            "Axis anisotropy: the kernel weight vanishes on this mesh at "
            "T = %.6g K, gamma = %.4e eV, so the aspect ratio cannot be "
            "estimated. Is mu inside the band range?", temperature, gamma,
        )
        return np.zeros(3)

    ncomp = moments.shape[2]
    d = np.zeros(3)
    for ax in range(3):
        if not dims[ax]:
            continue
        b = kvec[ax]
        q = (b[0] ** 2 * moments[:, :, 0]
             + b[1] ** 2 * moments[:, :, 1]
             + b[2] ** 2 * moments[:, :, 2])
        if ncomp >= 6:
            q = q + 2.0 * (b[0] * b[1] * moments[:, :, 3]
                           + b[0] * b[2] * moments[:, :, 4]
                           + b[1] * b[2] * moments[:, :, 5])
        d[ax] = np.sqrt(max(float(np.sum(omega * q)) / total, 0.0))

    if symop is not None:
        raw = d.copy()
        for orbit in _symmetry_axis_orbits(symop, dims):
            active = [a for a in orbit if dims[a]]
            if len(active) > 1:
                d[active] = float(np.mean(d[active]))
        # A large disagreement means the stored point group relates axes that
        # the BAND STRUCTURE does not, i.e. /.unitcell/symop overstates the
        # symmetry of the Hamiltonian.  ltb refuses to build an irreducible
        # mesh in that case, but a Wannier or DFT-derived file can still carry
        # it, and the irreducible reduction would then be invalid.
        scale = np.max(np.abs(raw)) or 1.0
        if np.max(np.abs(raw - d)) > 0.1 * scale:
            logger.warning(
                "The stored point group averages axes whose measured energy "
                "anisotropy differs substantially (%s before symmetrisation, "
                "%s after). The symmetry operations in /.unitcell/symop relate "
                "directions that the band structure does not, so this mesh's "
                "irreducible reduction may be invalid. Verify the symmetry, or "
                "work on a reducible grid.",
                np.array2string(raw, precision=4),
                np.array2string(d, precision=4),
            )
    return d


def suggest_mesh_ratio(
    d:        np.ndarray,
    dims:     np.ndarray,
    round_up: bool = True,
) -> "tuple[np.ndarray, np.ndarray]":
    """Raw and rounded-up axis ratios from :func:`axis_anisotropy`.

    Returns ``(raw, suggested)``, both normalised so the smallest ACTIVE axis
    is 1.  The suggestion rounds each ratio UP to the next integer, because
    the estimator is a lower bound (see :func:`axis_anisotropy`) and the
    optimum is broad, so over-shooting is the cheap direction to err in.

    Parity is deliberately NOT imposed here.  An odd division count misses the
    zone-boundary plane along that axis and is worth avoiding, but that is a
    property of the final divisions rather than of the ratio: forcing each
    ratio to be even and then renormalising halves it (a raw 6.97 came out as
    4:1 instead of 7:1).  Evenness is applied in :func:`recommend_density`,
    where actual division counts are formed.
    """

    d    = np.asarray(d, dtype=float).ravel()[:3]
    dims = np.asarray(dims).astype(bool).ravel()[:3]

    active = dims & (d > 0.0)
    raw = np.ones(3)
    sug = np.ones(3)
    if not np.any(active):
        return raw, sug

    base = float(np.min(d[active]))
    raw[active] = d[active] / base
    sug[active] = (np.maximum(np.ceil(raw[active] - 1e-9), 1.0) if round_up
                   else np.maximum(np.round(raw[active]), 1.0))
    sug[~dims] = 1.0
    return raw, sug


def recommend_density(
    nk:          "tuple[int, int, int]",
    error:       float,
    target:      float,
    dims:        np.ndarray,
    floor_error: "float | None" = None,
    multiple:    int = 1,
    max_growth:  float = 16.0,
) -> "tuple[tuple[int, int, int], bool, str]":
    """Suggest a uniform mesh that would reach *target* at the wide corner.

    The sum-rule metric on a UNIFORM mesh falls close to ``1/n`` (measured:
    graphene 0.1866 / 0.0708 / 0.0245 / 0.0098 for n = 48 / 96 / 192 / 300,
    cubic 0.2806 / 0.0326 / 0.0080 / 0.0050 for n = 8 / 16 / 24 / 32), so the
    factor by which every active axis must grow is estimated as
    ``error / target``.  That extrapolation is only used to make a
    suggestion; nothing depends on its accuracy.

    IMPORTANT -- the metric SATURATES.  Its floor is set by the band range and
    the trapezoid defect in the kernel tails, not by the mesh: about 0.0043
    for the cubic model and 0.002 for graphene.  A target at or below the
    floor is unreachable by ANY uniform mesh, and recommending an ever-larger
    nk would be a lie.  When *floor_error* is supplied and the target is not
    comfortably above it, this returns ``reachable = False`` and an
    explanation instead of a recommendation.

    Returns ``(nk, reachable, message)``.
    """

    nk   = tuple(int(v) for v in nk)
    dims = np.asarray(dims).astype(bool).ravel()[:3]

    if error <= target:
        return nk, True, "already at or below the target"

    if floor_error is not None and target <= floor_error * 1.2:
        return (
            nk, False,
            "the metric on this system saturates near %.4e -- set by the band "
            "range and the kernel tails, not by the mesh -- so the target "
            "%.4e cannot be reached by refining a uniform grid at all. Relax "
            "the target above roughly %.4e, or widen --energy_window."
            % (floor_error, target, floor_error * 1.2)
        )

    growth = min(float(error) / float(target), max_growth)
    out = []
    for ax in range(3):
        if not dims[ax] or nk[ax] <= 1:
            out.append(max(nk[ax], 1))
            continue
        n = int(np.ceil(nk[ax] * growth))
        if multiple > 1:
            n = int(np.ceil(n / multiple) * multiple)
        out.append(max(n, nk[ax] + 1))

    capped = ""
    if float(error) / float(target) > max_growth:
        capped = (" The required growth was capped at %.0fx per axis; expect "
                  "to repeat this step." % max_growth)
    return (tuple(out), True,
            "grow every active axis by about %.1fx.%s" % (growth, capped))


def infer_mesh_multiple(nk: "tuple[int, int, int]", dims: np.ndarray) -> int:
    """Largest small divisor shared by all active axes of *nk*.

    Preserving it keeps a recommended mesh commensurate with whatever the user
    already had -- which matters because high-symmetry points fall on the grid
    only for particular divisors.  Graphene needs n divisible by 3 to put K on
    the mesh; a parent without it stalls the refinement completely (measured:
    a 193x193 start froze at metric 1.339 and returned 9.79e3 against a
    4.43e5 reference).  Even counts are also preferred, since an odd count
    misses the zone-boundary plane along that axis.
    """

    dims = np.asarray(dims).astype(bool).ravel()[:3]
    active = [int(nk[a]) for a in range(3) if dims[a] and nk[a] > 1]
    if not active:
        return 1
    best = 1
    for m in (6, 4, 3, 2):
        if all(n % m == 0 for n in active):
            best = m
            break
    return best


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


def select_refinement_panels(
    band:               np.ndarray,
    temperature:        float,
    gamma:              float,
    chemical_potential: float,
    energy_window:      float = 0.1,
    defect_fraction:    float = 0.9,
    defect_rtol:        float = 1e-3,
) -> Tuple[np.ndarray, dict]:
    """Mark the band-axis panels that carry the bulk of the quadrature defect.

    This is the selection counterpart of :func:`compute_error`: the loop now
    optimises exactly the quantity it reports.  The previous criterion --
    ``detect_refinement_scale``, which flagged energies by the POINTWISE
    interpolation defect ``|analytic - interpolated|`` and turned their spread
    into a mu-centred radius -- measured something else, and did so with a
    systematic bias.  The kernel amplitude is ``1/(pi*gamma)`` at mu against
    ``gamma/(pi*eps^2)`` in the tails, a ratio of ``(eps/gamma)^2``, so a
    pointwise criterion always ranks the peak orders of magnitude above the
    tails no matter how coarse the tails are.  On graphene that pinned the
    window at its ``2*kernel_width`` floor from the fourth iteration onward
    while the residual sat between 0.01 and 0.1 eV, and the metric fell only
    as N^-0.26.

    Selecting whole PANELS rather than a radius also stops the refinement
    from being all-or-nothing in energy: a radius wide enough to reach a bad
    tail panel necessarily drags in every k-point closer to mu as well, which
    is what made the forced-wide-window experiment need 5.2e5 wedge points.

    Parameters
    ----------
    band : ndarray
        Sorted unique band energies (the band axis).
    temperature, gamma : float
        Kelvin and eV; the worst-case corner of the planned sweep.
    chemical_potential : float
        mu [eV]; used only for reporting and for the ``energy_window`` guard.
    energy_window : float
        Panels whose centre lies further than this from mu are never marked.
        Guards against chasing van Hove defects at the band edges, which carry
        real defect but no transport weight.
    defect_fraction : float
        Mark the largest-defect panels until their cumulative |d_i| reaches
        this fraction of the total defect inside the window.
    defect_rtol : float
        If the total defect is at or below this, report convergence by
        returning an all-False mask.

    Returns
    -------
    marked : ndarray (n_panels,) of bool
    info : dict
        Diagnostics for logging: ``n_marked``, ``n_panels``, ``total``,
        ``captured``, ``reach`` (largest |centre - mu| among marked panels).
    """

    mid, defects, _ = panel_defects(band, temperature, gamma)
    if defects.size == 0:
        return np.zeros(0, dtype=bool), {"n_marked": 0, "n_panels": 0,
                                         "total": 0.0, "captured": 0.0,
                                         "reach": 0.0}

    absd = np.abs(defects)
    # A panel counts as in-window if it OVERLAPS [mu - W, mu + W], not if its
    # centre happens to fall inside.  On a coarse mesh the panel straddling mu
    # can be far wider than the window -- graphene 24x24x1 has its first
    # energy above the Dirac point around 0.25 eV, so the single panel holding
    # the entire 81.9 defect has its centre at -0.125 eV and a centre test
    # discards it, reporting "converged" at an error of 81.9.
    band64 = np.asarray(band, dtype=np.float64).ravel()
    lo, hi = band64[:-1], band64[1:]
    mu     = float(chemical_potential)
    W      = float(energy_window)
    inside = (hi >= mu - W) & (lo <= mu + W)
    absd_in = np.where(inside, absd, 0.0)
    total = float(absd_in.sum())

    marked = np.zeros(defects.size, dtype=bool)
    info = {"n_marked": 0, "n_panels": int(defects.size),
            "total": total, "captured": 0.0, "reach": 0.0}

    if total <= float(defect_rtol):
        logger.info(
            "Panel defect inside the energy window is %.3e (<= defect_rtol=%.3e); "
            "no panel needs subdivision.", total, defect_rtol,
        )
        return marked, info

    order = np.argsort(absd_in)[::-1]
    cum   = np.cumsum(absd_in[order])
    need  = int(np.searchsorted(cum, float(defect_fraction) * total) + 1)
    need  = min(need, int(np.count_nonzero(absd_in)))
    take  = order[:need]
    marked[take] = True

    info["n_marked"] = int(need)
    info["captured"] = float(cum[need - 1]) if need > 0 else 0.0
    info["reach"]    = float(np.max(np.abs(mid[take] - float(chemical_potential))))
    return marked, info


def hotspot_mask_from_panels(
    energies: np.ndarray,
    band:     np.ndarray,
    marked:   np.ndarray,
) -> np.ndarray:
    """k-points that bracket at least one marked panel.

    A panel is bounded by two band-axis nodes, and every node is a band energy
    of at least one k-point.  Subdividing those k-points is what halves the
    panel width on the next iteration, so the hotspot set is exactly the
    k-points owning a node adjacent to a marked panel.
    """

    if marked.size == 0 or not marked.any():
        return np.zeros(energies.shape[0], dtype=bool)

    node_flag           = np.zeros(band.size, dtype=bool)
    node_flag[:-1]     |= marked
    node_flag[1:]      |= marked

    # every entry of *energies* is present in *band* by construction
    idx = np.searchsorted(np.asarray(band, dtype=np.float64),
                          np.asarray(energies, dtype=np.float64))
    idx = np.clip(idx, 0, band.size - 1)
    return np.any(node_flag[idx], axis=1)


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
    hotspot_mask:      "np.ndarray | None" = None,
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

    if hotspot_mask is not None:
        # Panel-defect selection (see select_refinement_panels): the hotspots
        # are the k-points bracketing the panels that carry the quadrature
        # error, which is not in general a mu-centred ball.
        mask    = np.asarray(hotspot_mask, dtype=bool)
        if mask.shape[0] != k_points.shape[0]:
            raise ValueError("hotspot_mask length {} does not match {} k-points"
                             .format(mask.shape[0], k_points.shape[0]))
        indices = np.where(mask)[0]
        if indices.size == 0:
            logger.warning("Hotspot mask selected no k-points.")
            return k_points.copy(), weights.copy(), cell_deltas.copy()
    else:
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
