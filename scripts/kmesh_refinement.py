"""
scripts/refine_kmesh.py

Generator-agnostic iterative k-mesh refinement loop.

This module contains no knowledge of tight-binding, Wannier, or any other
physical model.  It drives the refinement purely through the MeshGenerator
protocol defined in structure/meshgenerator.py.

Entry points for specific generators (ltb, lwann, ...) live in the
corresponding sub-packages under scripts/ and call run_refinement() after
constructing the appropriate generator object.
"""

from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import h5py
import numpy as np

from structure.meshgenerator import MeshGenerator
from structure.meshrefine import (
    BandData,
    MissingDependencyError,
    build_band_axis,
    compute_error,
    detect_refinement_scale,
    load_band_data,
    refine_kmesh,
    write_custom_mesh,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Parameters dataclass
# ---------------------------------------------------------------------------

@dataclass
class RefinementParams:
    """
    All tuning knobs for the refinement loop, collected in one place so that
    generator-specific entry points can build this object from their own
    argparse namespace without duplicating field names.

    Fields
    ------
    initial_hdf5 : Path
        HDF5 file produced by the first (coarse) generator run.
    chemical_potential : float
        mu used for hotspot detection and passed to the generator.
    gamma_min : float
        Minimum scattering rate (eV) for the df/da error metric.
    t_min : float
        Minimum temperature (K) for the df/da error metric.
    error_tol : float
        Stop when the df/da integration error falls below this value.
    max_iter : int
        Hard cap on refinement iterations.
    refinement_factor : int
        Number of subdivisions per k-axis for hotspot regions; should be odd.
    energy_window : float
        Half-width (eV) of the mu-centred window used by hotspot detection.
    workdir : Path
        Directory where intermediate mesh and HDF5 files are written.
    keep_intermediate : bool
        If True, do not delete intermediate files after each iteration.
    """
    initial_hdf5:        Path
    chemical_potential:  float
    gamma_min:           float
    t_min:               float
    error_tol:           float  = 5e-3
    max_iter:            int    = 10
    refinement_factor:   int    = 3
    energy_window:       float  = 0.1
    workdir:             Path   = Path('.')
    keep_intermediate:   bool   = False


# ---------------------------------------------------------------------------
# Shared refinement loop
# ---------------------------------------------------------------------------

def run_refinement(params: RefinementParams, generator: MeshGenerator) -> int:
    """
    Iteratively refine the k-mesh stored in *params.initial_hdf5*.

    On each iteration:
      1. Load band data from the current HDF5 file.
      2. Evaluate the df/da integration error.
      3. If the error is below *params.error_tol*, stop.
      4. Detect hotspot k-points near *params.chemical_potential*.
      5. Refine those k-points (subdivide) while preserving total weight.
      6. Call *generator.generate()* to produce the next HDF5 file.
      7. Optionally clean up intermediate files.

    Parameters
    ----------
    params : RefinementParams
        All loop parameters (tolerances, iteration cap, paths, …).
    generator : MeshGenerator
        Object that knows how to turn a (points, weights) mesh into an HDF5.

    Returns
    -------
    int
        0 on success, 1 on error (suitable as a process exit code).
    """
    assert isinstance(generator, MeshGenerator), (
        f"generator must implement the MeshGenerator protocol, got {type(generator)}"
    )

    workdir = params.workdir.resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    current_hdf5: Path = params.initial_hdf5.resolve()
    previous_output: Optional[Path] = None
    final_error: Optional[float] = None

    for iteration in range(params.max_iter):
        logger.info("--- Iteration %d ---", iteration)

        # ── 1. Load band data ─────────────────────────────────────────────
        with h5py.File(current_hdf5, 'r') as h5file:
            data: BandData = load_band_data(h5file)
            band_axis = build_band_axis(data.energies)

        # ── 2. Evaluate error ─────────────────────────────────────────────
        try:
            final_error = compute_error(band_axis, params.t_min, params.gamma_min)
        except MissingDependencyError as exc:
            logger.error("%s", exc)
            return 1

        logger.info("Iteration %d: error = %.6f", iteration, final_error)

        # ── 3. Convergence check ──────────────────────────────────────────
        if final_error <= params.error_tol:
            logger.info(
                "Target error reached (%.6f <= %.6f).", final_error, params.error_tol
            )
            break

        # ── 4. Detect hotspots ────────────────────────────────────────────
        try:
            hotspots, tolerance = detect_refinement_scale(
                band_axis,
                params.t_min,
                params.gamma_min,
                params.chemical_potential,
                energy_window=params.energy_window,
            )
        except MissingDependencyError as exc:
            logger.error("%s", exc)
            return 1

        if tolerance == 0.0:
            logger.warning("No significant hotspots detected; stopping early.")
            break

        logger.debug(
            "Hotspot count: %d, tolerance derived: %.6e", hotspots.size, tolerance
        )

        # ── 5. Refine mesh ────────────────────────────────────────────────
        refined_points, refined_weights = refine_kmesh(
            data,
            target_energy=params.chemical_potential,
            tolerance=tolerance,
            refinement_factor=params.refinement_factor,
        )

        if refined_points.shape[0] == data.k_points.shape[0]:
            logger.warning(
                "Refinement did not modify the mesh; stopping to avoid infinite loop."
            )
            break

        before_weight = np.sum(data.weights)
        after_weight  = np.sum(refined_weights)
        if not np.allclose(before_weight, after_weight, rtol=1e-10, atol=1e-12):
            logger.warning(
                "Weight conservation check failed: before=%s after=%s",
                before_weight, after_weight,
            )

        # ── 6. Write mesh file and call generator ─────────────────────────
        mesh_path   = workdir / f"custom_mesh_iter_{iteration + 1}.hdf5"
        output_path = workdir / f"refined_iter_{iteration + 1}.hdf5"

        if mesh_path.exists() and not params.keep_intermediate:
            mesh_path.unlink()

        write_custom_mesh(str(mesh_path), refined_points, refined_weights)

        try:
            generator.generate(refined_points, refined_weights, str(output_path))
        except Exception as exc:
            logger.error("Generator failed on iteration %d: %s", iteration, exc)
            return 1

        # ── 7. Clean up intermediate files ────────────────────────────────
        if not params.keep_intermediate:
            if mesh_path.exists():
                mesh_path.unlink()
            if (
                previous_output is not None
                and previous_output.exists()
                and previous_output != current_hdf5
            ):
                previous_output.unlink()

        previous_output = output_path
        current_hdf5    = output_path

    else:
        logger.warning(
            "Maximum iterations (%d) reached. Final error = %.6f",
            params.max_iter, final_error if final_error is not None else float('nan'),
        )

    logger.info("Refinement completed. Final mesh: %s", current_hdf5)
    return 0
