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


def trigamma_vec(z: np.ndarray) -> np.ndarray:
    """Vectorised trigamma using mpmath.polygamma with complex support."""

    _require_mpmath()
    z = np.asanyarray(z, dtype=np.complex128)
    out = np.empty_like(z, dtype=np.complex128)
    flat_in = z.ravel()
    flat_out = out.ravel()
    for idx, val in enumerate(flat_in):
        flat_out[idx] = complex(mp.polygamma(1, val))
    return out


def digamma_im(epsilon: np.ndarray, temperature: float, gamma: float) -> np.ndarray:
    """Imaginary part of the digamma-based auxiliary function."""

    beta = 1.0 / temperature
    z = 0.5 + (beta / (2.0 * np.pi)) * (gamma + 1j * epsilon)
    return 0.5 - (1.0 / np.pi) * np.imag(digamma(z))


def df_da(epsilon: np.ndarray, temperature: float, gamma: float) -> np.ndarray:
    """d/d(alpha) of the auxiliary function appearing in the error metric."""

    beta = 1.0 / temperature
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
    # by always writing the correct cell_deltas explicitly.
    cell_deltas:  np.ndarray   # (nkp, 3)


def load_band_data(h5file: "h5py.File") -> BandData:
    """Load band energies and k-mesh data from an HDF5 handle.

    cell_deltas is read from the ``/.kmesh/cell_deltas`` dataset when present.
    This dataset is written by :func:`write_custom_mesh` on every iteration
    from iteration 1 onward.

    For the initial coarse HDF5 (produced by ``ltb run`` / ``ltb refine``
    start, where cell_deltas is absent), the values are derived from the
    stored ``nkx`` / ``nky`` / ``nkz`` fields as ``1/nk_i``.  This path is
    only reached for the very first iteration.
    """

    energies = h5file["/energies"][:].astype(quad_dtype)
    if energies.ndim != 2 or energies.shape[1] < 1:
        raise ValueError("Energy dataset must be 2D with at least one band.")

    k_group      = h5file["/.kmesh"]
    k_points     = k_group["points"][:].astype(quad_dtype)
    weights      = k_group["weights"][:].astype(quad_dtype)
    multiplicity = k_group["multiplicity"][:].astype(quad_dtype)

    if "cell_deltas" in k_group:
        # Present from iteration 1 onward: written by write_custom_mesh.
        cell_deltas = k_group["cell_deltas"][:].astype(quad_dtype)
    else:
        # Initial coarse HDF5 from ltb run: nkx/nky/nkz are the true coarse
        # grid dimensions.  Derive uniform cell widths from them.
        nkx = int(k_group["nkx"][()])
        nky = int(k_group["nky"][()])
        nkz = int(k_group["nkz"][()])
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

    return BandData(
        energies=energies,
        k_points=k_points,
        weights=weights,
        multiplicity=multiplicity,
        cell_deltas=cell_deltas,
    )


def build_band_axis(bands: np.ndarray, max_bands: int = 2) -> np.ndarray:
    """Construct a monotonic energy axis from the lowest bands."""

    take = min(max_bands, bands.shape[1])
    subset = bands[:, :take]
    combined = np.sort(np.unique(subset.reshape(-1)))
    return combined.astype(quad_dtype)


def compute_error(band: np.ndarray, temperature: float, gamma: float) -> float:
    """Return |1 - int df_da d epsilon| for the supplied band grid."""

    dfda = df_da(band, temperature, gamma)
    integral = np.trapz(dfda, band)
    return float(np.abs(1.0 - integral))


def detect_refinement_scale(
    band: np.ndarray,
    temperature: float,
    gamma: float,
    chemical_potential: float,
    energy_window: float = 0.1,
    log_points: int = 110_000,
    uniform_points: int = 10_000,
) -> Tuple[np.ndarray, float]:
    """Identify undersampled energies and return their window-derived tolerance.

    Returns
    -------
    significant : ndarray
        Probe energies where the analytic df/da differs significantly from the
        linear interpolation of the coarse band grid.
    tolerance : float
        Energy window to use in :func:`refine_kmesh` for flagging hotspot
        k-points.  Equal to ``energy_window``, clipped to be at least 1e-6.
        This is the physically meaningful scale: the half-width of the region
        around ``chemical_potential`` that needs denser sampling.
    """

    dfda = df_da(band, temperature, gamma)

    band_min = float(np.min(band))
    band_max = float(np.max(band))

    uniform_band = np.linspace(band_min, band_max, num=uniform_points, endpoint=True, dtype=float)
    uniform_band = np.sort(np.unique(uniform_band))

    # Construct logarithmic sampling about the chemical potential; avoid zero
    # start by using max(|mu|, 1e-14) as the geomspace baseline.
    delta    = max(energy_window, 1e-8)
    baseline = max(abs(chemical_potential), 1e-14)
    geom     = np.geomspace(baseline, baseline + delta, num=log_points)
    log_band = np.sort(np.concatenate((chemical_potential + geom,
                                       chemical_potential - geom)))

    probe_band = np.sort(np.unique(np.concatenate((uniform_band, log_band))))

    interpolator = interp1d(band, dfda, kind="linear",
                            bounds_error=False, fill_value="extrapolate")
    interpolated = interpolator(probe_band)
    analytic     = df_da(probe_band, temperature, gamma)
    delta_arr    = np.abs(analytic - interpolated)

    if delta_arr.size == 0:
        return probe_band, 0.0

    threshold   = 0.5 * np.max(delta_arr)
    significant = probe_band[delta_arr >= threshold]

    if significant.size == 0:
        return probe_band, 0.0

    # The tolerance passed to refine_kmesh is the energy half-window around
    # chemical_potential within which k-points are considered hotspots.
    # energy_window is the physically correct scale for this: it is the
    # half-width of the df/da peak set by the user.  Using max|significant|*3
    # instead produced a window spanning the entire bandwidth for systems with
    # a Dirac point at mu (e.g. graphene), causing every k-point to be flagged.
    tolerance = max(float(energy_window), 1e-6)
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

    mask    = np.any(np.isclose(energies, target, atol=tol), axis=1)
    indices = np.where(mask)[0]

    if indices.size == 0:
        logger.warning(
            "No k-points within tolerance %.3e of target energy %.3e.",
            tolerance, target_energy,
        )
        return k_points.copy(), weights.copy(), cell_deltas.copy()

    keep_mask                 = np.ones(k_points.shape[0], dtype=bool)
    removed_weights           = []
    refined_points_collection = []
    refined_deltas_collection = []

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

    k_points_filtered = k_points[keep_mask]
    weights_filtered  = weights[keep_mask]
    deltas_filtered   = cell_deltas[keep_mask]

    total_removed = np.sum(removed_weights) if removed_weights else quad_dtype(0.0)

    if refined_points_collection:
        refined_points  = np.concatenate(refined_points_collection, axis=0)
        refined_deltas  = np.concatenate(refined_deltas_collection,  axis=0)
        n_refined       = refined_points.shape[0]
        refined_weights = np.full(n_refined,
                                  total_removed / quad_dtype(n_refined),
                                  dtype=quad_dtype)
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
