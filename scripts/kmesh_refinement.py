"""
scripts/kmesh_refinement.py

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
import re
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
    T_min : float
        Minimum temperature (K) for the df/da error metric.
    error_tol : float
        Stop when the df/da integration error falls below this value.
    max_iter : int
        Hard cap on refinement iterations.
    refinement_factor : int
        Number of subdivisions per k-axis for hotspot regions; should be odd.
    energy_window : float
        Half-width (eV) of the mu-centred window used by hotspot detection.
        Also used as the tolerance in refine_kmesh: k-points with any band
        energy within energy_window of mu are flagged as hotspots.
    workdir : Path
        Directory where intermediate mesh and HDF5 files are written.
    keep_intermediate : bool
        If True, do not delete intermediate files after each iteration.
    """
    initial_hdf5:        Path
    chemical_potential:  float
    gamma_min:           float
    T_min:               float
    error_tol:           float = 5e-3
    max_iter:            int   = 10
    refinement_factor:   int   = 3
    energy_window:       float = 0.1
    workdir:             Path  = Path('.')
    keep_intermediate:   bool  = False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_ITER_FILENAME_RE = re.compile(r"refined_iter_(\d+)\.hdf5$")


def _patch_refinement_metadata(
    output_path: Path,
    mesh_path: Path,
    iteration_index: int,
    parent: Path,
) -> None:
    """Copy /.kmesh/cell_deltas from *mesh_path* into *output_path* and stamp
    refinement provenance attributes.

    Background
    ----------
    ``write_custom_mesh`` writes ``cell_deltas`` to a lightweight mesh-only
    HDF5 (``custom_mesh_iter_N.hdf5``).  The generator then reads that file
    and produces a full output HDF5 (``refined_iter_N.hdf5``) via
    ``h5output``, which knows nothing about ``cell_deltas`` and therefore
    does not write it.

    The next iteration's ``load_band_data`` reads the generator output, finds
    no ``cell_deltas`` dataset, and would fall back to ``nkx/nky/nkz`` --
    which are always 1 for custom-mesh runs.  That would give zero cell
    widths, silently turning all subsequent subdivisions into no-ops
    (load_band_data now raises a hard error on that combination instead).

    The fix is to patch ``cell_deltas`` into the generator output immediately
    after it is written, before ``mesh_path`` is deleted.  Opening in append
    mode (``"a"``) leaves all other datasets untouched.

    Provenance attributes
    ---------------------
    ``refinement_iteration`` (int)
        Global 1-based index of the refinement iteration that produced this
        file.  Read back by :func:`_detect_refinement_state` so that a later
        ``ltb refine`` run starting from this file continues the numbering
        instead of restarting at 1 (and thereby colliding with -- or even
        overwriting -- earlier outputs, including the input file itself).
    ``refinement_parent`` (str)
        Path of the HDF5 file this iteration refined.  Informational.
    """
    with h5py.File(mesh_path,   "r") as src, \
         h5py.File(output_path, "a") as dst:
        cell_deltas = src["/.kmesh/cell_deltas"][:]
        kmesh = dst["/.kmesh"]
        if "cell_deltas" in kmesh:
            del kmesh["cell_deltas"]
        kmesh.create_dataset("cell_deltas", data=cell_deltas)
        dst.attrs["refinement_iteration"] = int(iteration_index)
        dst.attrs["refinement_parent"]    = str(parent)


def _detect_refinement_state(path: Path):
    """Classify *path* as fresh coarse mesh or (partially) refined mesh.

    Returns
    -------
    (is_continuation, last_iteration) : tuple[bool, Optional[int]]
        is_continuation : True if the file is an already-refined mesh, i.e.
            refinement should *continue* from it rather than start fresh.
        last_iteration : the 1-based index of the iteration that produced
            the file, if known (from the ``refinement_iteration`` attribute,
            or -- for files produced before that attribute existed -- parsed
            from a ``refined_iter_N.hdf5`` filename); None if unknown.

    Detection is hierarchical:
      1. ``refinement_iteration`` HDF5 attribute (robust against renames).
      2. Presence of ``/.kmesh/cell_deltas`` -> continuation; index possibly
         recovered from the filename.
      3. Otherwise: fresh coarse mesh.
    """
    with h5py.File(path, "r") as f:
        attr = f.attrs.get("refinement_iteration", None)
        has_deltas = "cell_deltas" in f["/.kmesh"]

    if attr is not None:
        return True, int(attr)

    if has_deltas:
        m = _ITER_FILENAME_RE.search(path.name)
        return True, (int(m.group(1)) if m else None)

    return False, None


def _max_existing_index(workdir: Path) -> int:
    """Highest N among refined_iter_N.hdf5 / custom_mesh_iter_N.hdf5 files
    already present in *workdir*; 0 if none.  Used so that a new run in a
    working directory containing outputs of earlier runs numbers its own
    outputs strictly after them, instead of interleaving into gaps.
    """
    pattern = re.compile(r"(?:refined|custom_mesh)_iter_(\d+)\.hdf5$")
    indices = [
        int(m.group(1))
        for f in workdir.iterdir()
        if (m := pattern.match(f.name))
    ]
    return max(indices, default=0)


def _reserve_iteration_paths(
    workdir: Path,
    index: int,
    forbidden: "set[Path]",
):
    """Return (mesh_path, output_path, index) for the next iteration such
    that neither file exists and neither resolves into *forbidden* (the
    input files of this run).  Existing files -- e.g. outputs of a previous
    refinement run in the same working directory -- are never overwritten;
    the index is bumped past them with a warning instead.
    """
    bumped_from = index
    while True:
        mesh_path   = workdir / f"custom_mesh_iter_{index}.hdf5"
        output_path = workdir / f"refined_iter_{index}.hdf5"
        collision = (
            mesh_path.exists()
            or output_path.exists()
            or mesh_path.resolve()   in forbidden
            or output_path.resolve() in forbidden
        )
        if not collision:
            break
        index += 1
    if index != bumped_from:
        logger.warning(
            "Output index bumped %d -> %d: refined_iter_%d.hdf5 / "
            "custom_mesh_iter_%d.hdf5 already exist in %s (or coincide with "
            "the input file) and will not be overwritten.",
            bumped_from, index, bumped_from, bumped_from, workdir,
        )
    return mesh_path, output_path, index


# ---------------------------------------------------------------------------
# Shared refinement loop
# ---------------------------------------------------------------------------

def run_refinement(params: RefinementParams, generator: MeshGenerator) -> int:
    """
    Iteratively refine the k-mesh stored in *params.initial_hdf5*.

    On each iteration:
      1. Load band data (energies, k-points, weights, cell_deltas) from the
         current HDF5 file.
      2. Evaluate the df/da integration error.
      3. If the error is below *params.error_tol*, stop.
      4. Detect hotspot k-points near *params.chemical_potential*.
      5. Refine those k-points (subdivide) while preserving total weight.
         Cell widths are taken from the stored cell_deltas array, so the
         subdivision is correct regardless of whether the starting mesh is
         reducible or irreducible and regardless of which iteration created
         each point.
      6. Call *generator.generate()* to produce the next HDF5 file.
      7. Patch cell_deltas and provenance metadata from the mesh file into
         the generator output so that they are available on the next
         iteration or a later continuation run (_patch_refinement_metadata).
      8. Optionally clean up intermediate files (only files created by this
         run are ever deleted; the input file is never touched).

    Continuation
    ------------
    *params.initial_hdf5* may itself be the output of a previous refinement
    run (e.g. ``refined_iter_3.hdf5``).  This is detected automatically from
    the ``refinement_iteration`` attribute / the ``cell_deltas`` dataset:
    per-point cell widths are then read from the file (NOT re-derived from
    coarse grid steps), and output numbering resumes after the input's
    iteration index.  Output files never overwrite existing files: if a name
    is taken, the index is bumped past it with a warning.

    Parameters
    ----------
    params : RefinementParams
    generator : MeshGenerator

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

    current_hdf5: Path            = params.initial_hdf5.resolve()
    final_error:  Optional[float] = None

    # files this run created and may therefore delete; the input file is
    # never in this set and is consequently never deleted or overwritten
    created_by_run: set = set()

    # ── Continuation detection ────────────────────────────────────────────
    try:
        is_continuation, last_iter = _detect_refinement_state(current_hdf5)
    except Exception as exc:
        logger.error("Cannot inspect %s: %s", current_hdf5, exc)
        return 1
    # resume after the input's iteration index, and never number below
    # refinement files already present in the working directory
    next_index = max(
        (last_iter + 1) if last_iter is not None else 1,
        _max_existing_index(workdir) + 1,
    )

    for iteration in range(params.max_iter):
        logger.info("--- Iteration %d ---", iteration)

        # ── 1. Load band data ─────────────────────────────────────────────
        try:
            with h5py.File(current_hdf5, 'r') as h5file:
                data: BandData = load_band_data(h5file)
                band_axis = build_band_axis(data.energies)

                if iteration == 0:
                    if is_continuation:
                        assert data.cell_deltas_source == "file", (
                            "continuation input must provide /.kmesh/cell_deltas"
                        )
                        logger.info(
                            "Initial mesh is an already-refined custom mesh "
                            "(%d k-points; produced by refinement iteration %s). "
                            "Continuing refinement: per-point cell widths read "
                            "from /.kmesh/cell_deltas, output numbering resumes "
                            "at %d.",
                            data.k_points.shape[0],
                            str(last_iter) if last_iter is not None else "unknown",
                            next_index,
                        )
                    else:
                        irreducible = bool(h5file['/.kmesh/irreducible'][()])
                        mesh_type   = "irreducible" if irreducible else "reducible"
                        logger.info(
                            "Initial mesh is %s (%d k-points). "
                            "Cell widths initialised from coarse grid steps 1/nk_i.",
                            mesh_type, data.k_points.shape[0],
                        )
        except ValueError as exc:
            # e.g. custom-mesh file without cell_deltas (see load_band_data)
            logger.error("%s", exc)
            return 1

        # ── 2. Evaluate error ─────────────────────────────────────────────
        try:
            final_error = compute_error(band_axis, params.T_min, params.gamma_min)
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
                params.T_min,
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

        logger.info(
            "Hotspot tolerance: %.4e eV  (energy_window=%.4e, "
            "%d probe energies flagged as undersampled).",
            tolerance, params.energy_window, hotspots.size,
        )

        # ── 5. Refine mesh ────────────────────────────────────────────────
        refined_points, refined_weights, refined_deltas = refine_kmesh(
            data,
            target_energy=params.chemical_potential,
            tolerance=tolerance,
            refinement_factor=params.refinement_factor,
        )

        n_before = data.k_points.shape[0]
        n_after  = refined_points.shape[0]

        if n_after == n_before:
            logger.warning(
                "Refinement did not modify the mesh; stopping to avoid "
                "infinite loop."
            )
            break

        logger.info("Mesh size: %d → %d k-points.", n_before, n_after)

        before_weight = np.sum(data.weights)
        after_weight  = np.sum(refined_weights)
        if not np.allclose(before_weight, after_weight, rtol=1e-10, atol=1e-12):
            logger.warning(
                "Weight conservation check failed: before=%s after=%s",
                before_weight, after_weight,
            )

        # ── 6. Write mesh file and call generator ─────────────────────────
        # Collision-free naming: never overwrite existing files (outputs of
        # earlier runs, or the input file of a continuation run in the same
        # working directory).  next_index resumes after the input's iteration
        # for continuation runs.
        mesh_path, output_path, next_index = _reserve_iteration_paths(
            workdir, next_index, forbidden={params.initial_hdf5.resolve()},
        )

        write_custom_mesh(str(mesh_path), refined_points, refined_weights,
                          refined_deltas)
        created_by_run.add(mesh_path)

        try:
            generator.generate(refined_points, refined_weights, str(output_path))
        except Exception as exc:
            logger.error("Generator failed on iteration %d: %s", iteration, exc)
            return 1
        created_by_run.add(output_path)

        # ── 7. Patch cell_deltas and provenance into the generator output ─
        # generator.generate() calls h5output which does not know about
        # cell_deltas.  Copy it from mesh_path into output_path so that the
        # next iteration's (or a later continuation run's) load_band_data
        # finds it directly and never falls back to the nkx/nky/nkz path
        # (which always gives 1 for custom meshes).  Also stamp the global
        # iteration index so continuation runs resume the numbering.
        try:
            _patch_refinement_metadata(output_path, mesh_path,
                                       iteration_index=next_index,
                                       parent=current_hdf5)
        except Exception as exc:
            logger.error("Failed to patch refinement metadata into %s: %s",
                         output_path, exc)
            return 1

        # ── 8. Clean up intermediate files ────────────────────────────────
        # Only files created by this run are ever deleted; the input file
        # (params.initial_hdf5) is never in created_by_run.
        if not params.keep_intermediate:
            if mesh_path.exists():
                mesh_path.unlink()
            if current_hdf5 in created_by_run and current_hdf5.exists():
                # output of the previous iteration of this run, now superseded
                current_hdf5.unlink()

        current_hdf5 = output_path
        next_index  += 1

    else:
        logger.warning(
            "Maximum iterations (%d) reached. Final error = %.6f",
            params.max_iter,
            final_error if final_error is not None else float('nan'),
        )

    logger.info("Refinement completed. Final mesh: %s", current_hdf5)
    return 0
