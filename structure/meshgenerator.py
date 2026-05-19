"""
structure/meshgenerator.py

Defines the MeshGenerator protocol: the interface that any k-mesh-aware
HDF5 generator must satisfy to work with the shared refinement loop in
scripts/refine_kmesh.py.

Adding a new generator (e.g. for lwann or lhk) requires only implementing
a class with a generate() method matching the signature below -- no changes
to the refinement loop or to existing generators are needed.
"""

from __future__ import annotations

from typing import runtime_checkable

import numpy as np

try:
    from typing import Protocol          # Python >= 3.8
except ImportError:                      # pragma: no cover
    from typing_extensions import Protocol


@runtime_checkable
class MeshGenerator(Protocol):
    """
    Anything that can evaluate a physical model on a custom k-mesh and write
    the result to an HDF5 file.

    Implementations
    ---------------
    structure.generators.ltb_gen.LtbGenerator   -- tight-binding via ltb
    (future) structure.generators.lwann_gen.WannGenerator
    (future) structure.generators.lhk_gen.HkGenerator
    """

    def generate(
        self,
        mesh_points: np.ndarray,
        mesh_weights: np.ndarray,
        output_path: str,
    ) -> None:
        """
        Evaluate the model on the supplied k-mesh and write the output HDF5.

        Parameters
        ----------
        mesh_points : ndarray, shape (N, 3)
            k-points in fractional coordinates.
        mesh_weights : ndarray, shape (N,)
            Integration weights.  Weight conservation across refinement
            iterations is the caller's responsibility.
        output_path : str
            Destination HDF5 file path.  Existing files are overwritten.
        """
        ...
