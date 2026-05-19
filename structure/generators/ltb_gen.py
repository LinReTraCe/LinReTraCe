"""
structure/generators/ltb_gen.py

MeshGenerator implementation for tight-binding (ltb) calculations.

Wraps the TightBinding class so that the iterative mesh-refinement loop in
scripts/refine_kmesh.py can drive it without knowing any TB-specific details.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import numpy as np

from structure.inout import h5output
from structure.tb import TightBinding

logger = logging.getLogger(__name__)


class LtbGenerator:
    """
    Generate a LinReTraCe-ready HDF5 file from a tight-binding model evaluated
    on a caller-supplied custom k-mesh.

    Parameters
    ----------
    tb_file : str or Path
        Tight-binding input file (Wannier-style, with 'begin hopping' etc.).
    filling : float
        Number of electrons in the system.
    mu : float or None
        Fixed chemical potential.  If None the Fermi level is determined
        automatically from *filling*.
    mushift : bool
        Shift all energies so that mu = 0 after diagonalisation.
    corronly : bool
        Use only the multi-orbital Peierls correction term (set derivative
        term to zero).
    vector : bool
        Also save the full Hamiltonian and transformation matrices to the
        output HDF5.
    intra : float or None
        Override all intra-band optical matrix elements with this value.
    inter : float or None
        Override all inter-band optical matrix elements with this value.
    intraonly : bool
        Discard inter-band elements before writing (sets opticfull = False).
    """

    def __init__(
        self,
        tb_file: str | Path,
        filling: float = 2.0,
        mu: float | None = None,
        mushift: bool = False,
        corronly: bool = False,
        vector: bool = False,
        intra: float | None = None,
        inter: float | None = None,
        intraonly: bool = False,
    ) -> None:
        self.tb_file   = str(tb_file)
        self.filling   = filling
        self.mu        = mu
        self.mushift   = mushift
        self.corronly  = corronly
        self.vector    = vector
        self.intra     = intra
        self.inter     = inter
        self.intraonly = intraonly

    # ------------------------------------------------------------------
    # Public interface expected by the MeshGenerator protocol
    # ------------------------------------------------------------------

    def generate(
        self,
        mesh_points: np.ndarray,
        mesh_weights: np.ndarray,
        output_path: str | Path,
    ) -> None:
        """
        Evaluate the tight-binding model on *mesh_points* / *mesh_weights*
        and write the result to *output_path*.

        Parameters
        ----------
        mesh_points : ndarray, shape (N, 3)
            k-points in fractional coordinates.
        mesh_weights : ndarray, shape (N,)
            Integration weights (must sum to the same total as the original
            mesh — weight conservation is the caller's responsibility).
        output_path : str or Path
            Destination HDF5 file.  Any existing file at this path is
            overwritten.
        """
        output_path = str(output_path)

        # Validate mesh arrays before touching TightBinding
        mesh_points  = np.asarray(mesh_points,  dtype=float)
        mesh_weights = np.asarray(mesh_weights, dtype=float).ravel()
        if mesh_points.ndim != 2 or mesh_points.shape[1] != 3:
            raise ValueError(
                "mesh_points must be an (N, 3) array of fractional coordinates."
            )
        if mesh_weights.shape[0] != mesh_points.shape[0]:
            raise ValueError(
                "mesh_weights must have the same length as mesh_points."
            )

        logger.info(
            "LtbGenerator: evaluating TB model on %d k-points -> %s",
            mesh_points.shape[0], output_path,
        )

        # Build TightBinding object (nkx/nky/nkz are placeholders for custom meshes)
        tb = TightBinding(nkx=1, nky=1, nkz=1, irreducible=False, kshift=False)

        tb.setCustomKmesh(mesh_points, mesh_weights)
        logger.info("Custom k-mesh set: %d k-points.", tb.nkp)

        tb.computeData(
            tbfile   = self.tb_file,
            charge   = self.filling,
            mu       = self.mu,
            mushift  = self.mushift,
            corronly = self.corronly,
            vector   = self.vector,
        )

        # Optional overrides of optical matrix elements
        if self.intra is not None:
            tb.setDiagonal(self.intra)
            logger.info("Intra-band elements set to %.6f.", self.intra)

        if self.inter is not None:
            tb.setOffDiagonal(self.inter)
            logger.info("Inter-band elements set to %.6f.", self.inter)

        if self.intraonly:
            tb.bopticfull = False
            tb.opticfull  = False
            logger.info("Inter-band elements will not be written (intraonly).")

        h5output(output_path, tb, tb)
        size_mb = os.stat(output_path).st_size / 1e6
        logger.info("Written %s (%.2f MB).", output_path, size_mb)
