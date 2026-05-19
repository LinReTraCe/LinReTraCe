"""Utilities for iterative k-mesh refinement based on energy-space error estimates."""

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
    energies: np.ndarray
    k_points: np.ndarray
    weights: np.ndarray
    multiplicity: np.ndarray


def load_band_data(h5file: "h5py.File") -> BandData:
    """Load band energies and k-mesh data from an HDF5 handle."""

    energies = h5file["/energies"][:].astype(quad_dtype)
    if energies.ndim != 2 or energies.shape[1] < 1:
        raise ValueError("Energy dataset must be 2D with at least one band.")

    k_group = h5file["/.kmesh"]
    k_points = k_group["points"][:].astype(quad_dtype)
    weights = k_group["weights"][:].astype(quad_dtype)
    multiplicity = k_group["multiplicity"][:].astype(quad_dtype)

    return BandData(energies=energies, k_points=k_points, weights=weights, multiplicity=multiplicity)


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
    """Identify undersampled energies and return their window-derived tolerance."""

    dfda = df_da(band, temperature, gamma)

    band_min = float(np.min(band))
    band_max = float(np.max(band))

    uniform_band = np.linspace(band_min, band_max, num=uniform_points, endpoint=True, dtype=float)
    uniform_band = np.sort(np.unique(uniform_band))

    # Construct logarithmic sampling about the chemical potential; avoid zero start.
    delta = max(energy_window, 1e-8)
    baseline = max(abs(chemical_potential), 1e-14)
    geom = np.geomspace(baseline, baseline + delta, num=log_points)
    log_band = np.sort(np.concatenate((chemical_potential + geom, chemical_potential - geom)))

    probe_band = np.sort(np.unique(np.concatenate((uniform_band, log_band))))

    interpolator = interp1d(band, dfda, kind="linear", bounds_error=False, fill_value="extrapolate")
    interpolated = interpolator(probe_band)
    analytic = df_da(probe_band, temperature, gamma)
    delta_arr = np.abs(analytic - interpolated)

    if delta_arr.size == 0:
        return probe_band, 0.0

    threshold = 0.5 * np.max(delta_arr)
    significant = probe_band[delta_arr >= threshold]

    if significant.size == 0:
        return probe_band, 0.0

    refine_window = np.max(np.abs(significant)) * 3.0
    refine_window = max(float(refine_window), 1e-6)
    return significant.astype(quad_dtype), refine_window


def _merge_duplicates_with_tolerance(
    points: np.ndarray,
    weights: np.ndarray,
    tol: float,
) -> Tuple[np.ndarray, np.ndarray]:
    if points.size == 0:
        return points, weights

    scaled = np.round(points / tol) * tol
    dtype = scaled.dtype
    structured = scaled.view([(f"f{i}", dtype) for i in range(scaled.shape[1])])
    _, idx, inv = np.unique(structured, return_index=True, return_inverse=True)
    merged_points = points[idx]
    merged_weights = np.zeros(len(idx), dtype=weights.dtype)
    np.add.at(merged_weights, inv, weights)
    return merged_points, merged_weights


def refine_kmesh(
    band_data: BandData,
    target_energy: float,
    tolerance: float,
    refinement_factor: int = 3,
    uniqueness_tol: float = 1e-12,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return refined (points, weights) preserving total weight.
    
    The refinement uses a cell-centred subdivision: each hotspot k-point
    represents a cell of linear size d (its nearest-neighbour spacing).
    That cell is split into refinement_factor^3 sub-cells, with new points
    placed at sub-cell centres:
 
        k_new = k_parent + (i - (n-1)/2) * d/n,   i = 0 ... n-1
 
    For odd n one sub-cell centre coincides exactly with the parent point,
    so the coarse mesh is always contained in the refined mesh and the
    subdivision is symmetric around the parent.  Even refinement_factor
    values are therefore not allowed and are incremented by 1 with a warning.
    """

    # Enforce odd refinement_factor
    if refinement_factor % 2 == 0:
        refinement_factor += 1
        logger.warning(
            "refinement_factor must be odd for symmetric cell-centred subdivision. "
            "Incrementing to %d.", refinement_factor
        )
    
    energies = band_data.energies
    k_points = band_data.k_points
    weights = band_data.weights

    tol = quad_dtype(tolerance)
    target = quad_dtype(target_energy)

    mask = np.any(np.isclose(energies, target, atol=tol), axis=1)
    indices = np.where(mask)[0]

    if indices.size == 0:
        logger.warning("No k-points within tolerance %.3e of target energy %.3e.", tolerance, target_energy)
        return k_points.copy(), weights.copy()

    keep_mask = np.ones(k_points.shape[0], dtype=bool)
    removed_weights = []
    refined_points_collection = []

    one = quad_dtype(1.0)

    for idx in indices:
        keep_mask[idx] = False
        removed_weights.append(weights[idx])

        k_pt = k_points[idx]
        diffs = []
        for axis in range(3):
            axis_diff = np.abs(k_points[:, axis] - k_pt[axis])
            axis_diff = axis_diff[axis_diff > quad_dtype(1e-33)]
            diffs.append(axis_diff.min() if axis_diff.size > 0 else quad_dtype(0.0))

        deltas = np.array(diffs, dtype=quad_dtype)
        step = deltas / quad_dtype(refinement_factor)

        # Cell-centred subdivision: place n points symmetrically around the
        # parent at offsets (i - (n-1)/2) * d/n for i = 0 ... n-1.
        # For odd n the central offset is 0, so the parent is always included.
        n = refinement_factor
        offsets = (np.arange(n, dtype=quad_dtype) - quad_dtype(n - 1) / 2)
        refined_axes = [
            k_pt[axis] + offsets * (deltas[axis] / quad_dtype(n))
            if deltas[axis] > 0 else np.array([k_pt[axis]], dtype=quad_dtype)
            for axis in range(3)
        ]

        # old scheme
        #refined_axes = [
        #    np.linspace(k_pt[axis] - deltas[axis] / 2.0, k_pt[axis] + deltas[axis] / 2.0,
        #                refinement_factor, endpoint=False, dtype=float)
        #    if deltas[axis] > 0 else np.array([float(k_pt[axis])])
        #    for axis in range(3)
        #]

        mesh = np.meshgrid(*refined_axes, indexing="ij")
        refined = np.column_stack([axis.ravel() for axis in mesh]).astype(quad_dtype)
        refined %= one
        refined[refined >= (one - quad_dtype(1e-15))] = quad_dtype(0.0)
        refined_points_collection.append(refined)

    k_points_filtered = k_points[keep_mask]
    weights_filtered = weights[keep_mask]

    total_removed = np.sum(removed_weights) if removed_weights else quad_dtype(0.0)
    refined_points = np.concatenate(refined_points_collection, axis=0) if refined_points_collection else np.empty((0, 3), dtype=quad_dtype)

    if refined_points.size:
        n_refined = refined_points.shape[0]
        refined_weights = np.full(n_refined, total_removed / quad_dtype(n_refined), dtype=quad_dtype)
    else:
        refined_weights = np.empty((0,), dtype=quad_dtype)

    merged_points = np.concatenate((k_points_filtered, refined_points), axis=0)
    merged_weights = np.concatenate((weights_filtered, refined_weights), axis=0)

    merged_points, merged_weights = _merge_duplicates_with_tolerance(
        merged_points, merged_weights, tol=quad_dtype(uniqueness_tol)
    )

    return merged_points, merged_weights


def write_custom_mesh(path: str, points: np.ndarray, weights: np.ndarray) -> None:
    import h5py

    with h5py.File(path, "w") as h5file:
        km_group = h5file.create_group(".kmesh")
        km_group.create_dataset("points", data=points.astype(float))
        km_group.create_dataset("weights", data=weights.astype(float))
