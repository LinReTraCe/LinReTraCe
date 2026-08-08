"""
structure/generators/lwann_gen.py

MeshGenerator implementation for Wannier90 (lwann) calculations.

Wraps Wannier90Calculation so that the iterative mesh-refinement loop in
scripts/kmesh_refinement.py can drive it without knowing any Wannier-specific
details.  The counterpart for tight-binding models is
structure/generators/ltb_gen.py.

Baseline implementation
-----------------------
Every call to :meth:`WannGenerator.generate` builds a fresh
Wannier90Calculation, re-reads the Wannier90 folder, and evaluates the model
on the *entire* mesh -- exactly like LtbGenerator does for tight-binding.
This is deliberately the simplest thing that can work, so that the timings
logged below measure the real cost split (read / compute / write) and can be
used to decide whether an incremental generator (evaluating only the newly
created child k-points and merging with the parent's stored arrays) is worth
the extra protocol surface.
"""

from __future__ import annotations

import logging
import os
import time
from pathlib import Path

import numpy as np

from structure.wannier import Wannier90Calculation

logger = logging.getLogger(__name__)


class WannGenerator:
    """
    Generate a LinReTraCe-ready HDF5 file from a Wannier90 calculation
    evaluated on a caller-supplied custom k-mesh.

    Parameters
    ----------
    directory : str or Path
        Folder of the Wannier90 calculation (as passed to ``lwann``).
    charge : float or None
        Number of electrons in the projected bands.  If None, the charge is
        determined at *mu* and rounded to the nearest integer, mirroring the
        behaviour of ``lwann`` without ``--charge``.
    mu : float or None
        Chemical potential.  Used when *charge* is None.  When *charge* is
        given, mu is instead determined from the charge on each mesh and may
        therefore drift slightly between refinement iterations.
    soc : bool
        Wannier90 calculation includes spin-orbit coupling (only needed for
        spin-unpolarised calculations).
    peierlscorrection : bool
        Include the multi-atomic (generalised) Peierls correction term.
    intraonly : bool
        Discard inter-band elements before writing.  Also passed to
        computeHamiltonian, which then never allocates the full
        (nband x nband) moment arrays -- the dominant memory term.
    kblock : int or None
        k-points evaluated per pass inside computeHamiltonian.  None sizes
        the block from *memory_budget_gb*.
    memory_budget_gb : float or None
        Budget [GB] for the transient per-block arrays.  None uses a quarter
        of the detected available memory.
    dims : sequence of 3 bool or None
        Dimensionality of the parent calculation (.unitcell/dims of the coarse
        HDF5).  Adopted verbatim so that the refined file is treated exactly
        like the mesh it was derived from -- this matters for setups whose
        k-mesh is three-dimensional while the Wannier projection is not (e.g.
        a single planar d_x2-y2 orbital in a layered cuprate), where inferring
        the dimensions from the k-points would give the wrong answer.  If
        None, the dimensions are inferred from the k-point spread.
    """

    def __init__(
        self,
        directory,
        charge: "float | None" = None,
        mu: "float | None" = None,
        soc: bool = False,
        peierlscorrection: bool = True,
        intraonly: bool = False,
        dims=None,
        kblock: "int | None" = None,
        memory_budget_gb: "float | None" = None,
    ):
        self.directory         = Path(directory)
        self.charge            = charge
        self.mu                = mu
        self.soc               = bool(soc)
        self.peierlscorrection = bool(peierlscorrection)
        self.intraonly         = bool(intraonly)
        self.dims              = dims
        self.kblock            = kblock
        self.memory_budget_gb  = memory_budget_gb

        if self.charge is None and self.mu is None:
            raise ValueError("WannGenerator requires either charge or mu.")

        if not self.directory.is_dir():
            raise IOError(
                "WannGenerator needs a Wannier90 folder, got: {}".format(self.directory)
            )

    def generate(
        self,
        mesh_points:  np.ndarray,
        mesh_weights: np.ndarray,
        output_path:  str,
    ) -> None:
        """Evaluate the Wannier model on *mesh_points* and write the HDF5."""

        output_path = str(output_path)

        # The refinement loop works in extended precision; the Wannier code
        # and h5output both require float64.
        mesh_points  = np.asarray(mesh_points,  dtype=np.float64)
        mesh_weights = np.asarray(mesh_weights, dtype=np.float64).ravel()
        if mesh_points.ndim != 2 or mesh_points.shape[1] != 3:
            raise ValueError(
                "mesh_points must be an (N, 3) array of fractional coordinates."
            )
        if mesh_weights.shape[0] != mesh_points.shape[0]:
            raise ValueError(
                "mesh_weights must have the same length as mesh_points."
            )

        logger.info(
            "WannGenerator: evaluating Wannier90 model on %d k-points -> %s",
            mesh_points.shape[0], output_path,
        )

        t_read0 = time.perf_counter()
        ham = Wannier90Calculation(str(self.directory), self.charge, self.soc)
        ham.readData()
        t_read = time.perf_counter() - t_read0

        ham.setCustomKmesh(mesh_points, mesh_weights, dims=self.dims)

        t_comp0 = time.perf_counter()
        ham.computeHamiltonian(
            peierlscorrection = self.peierlscorrection,
            kblock            = self.kblock,
            memory_budget_gb  = self.memory_budget_gb,
            intraonly         = self.intraonly,
        )
        t_comp = time.perf_counter() - t_comp0

        t_out0 = time.perf_counter()
        if self.charge is not None:
            # charge fixed -> mu follows from the charge on this mesh
            ham.outputData(output_path, intraonly=self.intraonly)
        else:
            # mu fixed -> charge is determined at mu and rounded (as in lwann)
            ham.outputData(output_path, mu=self.mu, intraonly=self.intraonly)
        t_out = time.perf_counter() - t_out0

        size_mb = os.stat(output_path).st_size / 1e6
        logger.info(
            "Written %s (%.2f MB). Timings [s]: read %.2f | compute %.2f | "
            "output %.2f (%d k-points, %d projections).",
            output_path, size_mb, t_read, t_comp, t_out,
            ham.nkp, ham.nproj,
        )
