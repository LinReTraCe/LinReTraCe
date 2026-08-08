#!/usr/bin/env python3
"""
scripts/lwann/refine_kmesh.py

lwann entry point for iterative k-mesh refinement.

Parses command-line arguments, constructs a WannGenerator, translates the
arguments into a RefinementParams object, and hands both to the shared
generator-agnostic loop in scripts/kmesh_refinement.py.  The loop itself is
identical to the one driven by `ltb refine`; only the generator differs.

Usage
-----
    lwann refine initial.hdf5 wannier_dir mu gamma_min T_min [options]
    python -m scripts.lwann.refine_kmesh initial.hdf5 wannier_dir mu
                                          gamma_min T_min [options]

Positional arguments
--------------------
    initial_hdf5        Initial Wannier90 HDF5 file (must end with .hdf5).
                        May be either a coarse mesh from a plain `lwann` run
                        OR the output of a previous refinement run
                        (refined_iter_N.hdf5): the latter is detected
                        automatically and refinement then CONTINUES from it --
                        per-point cell widths are read from the stored
                        /.kmesh/cell_deltas and output numbering resumes after
                        N.  Existing files in the working directory are never
                        overwritten.
    wannier_dir         Folder of the Wannier90 calculation.
    chemical_potential  Chemical potential (mu) used for df/da and, unless
                        --charge is given, for the generator.
    gamma_min           Minimum gamma parameter [eV] for df/da.
    T_min               Minimum temperature [K] for df/da.

Optional arguments
------------------
    --error_tol ERR         Target df/da error threshold (default: 0.005).
    --plateau_tol PERCENT   Error-plateau detection threshold (default: 5).
    --max_iter N            Maximum refinement iterations (default: 10).
    --refinement_factor F   Subdivisions per axis for refined hotspots (default: 3).
    --energy_window E       mu-centred energy window ceiling (default: 0.1).
    --charge C              Number of electrons in the projected bands.
    --soc                   Wannier90 calculation includes spin-orbit coupling.
    --nocorrection          Skip the Peierls correction term.
    --intraonly             Discard inter-band elements before writing.
    --kblock N              k-points per pass when building H(k).
    --max-memory GB         Budget for transient per-block arrays.
    --keep_intermediate     Keep intermediate mesh and HDF5 files.
    --workdir DIR           Directory for intermediate and output files (default: cwd).
    --openmp [N]            Number of OpenMP/BLAS threads; omit N for all cores.
    --verbose               Enable debug logging.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure the package root is on sys.path when this script is run directly
# by absolute path. scripts/lwann/ sits two levels below the root.
_root = Path(__file__).resolve().parents[2]
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

# ── OpenMP / BLAS thread count ────────────────────────────────────────────────
# MUST happen before numpy (or any BLAS-linked library) is imported.
# NOTE: openmp_utils currently lives under scripts/ltb/ but is entirely
# generator-agnostic; it should move to scripts/ in the planned refactor.
from scripts.ltb.openmp_utils import preparse_openmp, add_openmp_argument, openmp_info_line
_openmp_ncores = preparse_openmp()
# ─────────────────────────────────────────────────────────────────────────────

import argparse
import logging

from structure.generators.lwann_gen import WannGenerator
from scripts.kmesh_refinement import (RefinementParams, read_source_dims,
                                       read_source_symmetry, run_refinement)

LOG_FORMAT = "%(levelname)s: %(message)s"
logger = logging.getLogger(__name__)

# See scripts/ltb/refine_kmesh.py for the rationale.
T_MIN_FLOOR_K = 5e-4


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog='lwann refine',
        description='Iteratively refine a LinReTraCe k-mesh from a Wannier90 '
                    'calculation towards a target error threshold.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Usage via wrapper (recommended)
--------------------------------
  lwann refine initial.hdf5 wannier_dir mu gamma_min T_min [options]

Usage via direct call (expert)
--------------------------------
  python -m scripts.lwann.refine_kmesh initial.hdf5 wannier_dir mu gamma_min T_min [options]

Positional arguments
--------------------
  initial_hdf5        Initial Wannier90 HDF5 file (must end with .hdf5).
                      A previous refinement output (refined_iter_N.hdf5) is
                      detected automatically -> refinement continues from it.
  wannier_dir         Folder of the Wannier90 calculation.
  chemical_potential  Chemical potential (mu) used for df/da and the generator.
  gamma_min           Minimum gamma parameter for df/da.
  T_min               Minimum temperature for df/da.
""",
    )

    # Positional arguments
    parser.add_argument("initial_hdf5", type=Path,
                        help="Initial Wannier90 HDF5 file (must end with .hdf5). "
                             "May be a coarse mesh from a plain 'lwann' run or the "
                             "output of a previous refinement run; the latter is "
                             "detected automatically and refinement continues from it.")
    parser.add_argument("wannier_dir", type=Path,
                        help="Folder of the Wannier90 calculation.")
    parser.add_argument("chemical_potential", type=float,
                        help="Chemical potential [eV] used for refinement and, "
                             "unless --charge is given, for the generator.")
    parser.add_argument("gamma_min", type=float,
                        help="Minimum gamma parameter [eV] for df/da.")
    parser.add_argument("T_min", type=float,
                        help="Minimum temperature [K] for df/da; enters as "
                        "beta = 1/(kB*T) in eV^-1. If 0 is supplied it is "
                        "raised to %g K, since beta diverges at T=0."
                        % T_MIN_FLOOR_K)

    # Refinement loop options (identical to ltb refine)
    parser.add_argument("--error_tol", type=float, default=5e-3,
                        help="Target error threshold. Default: 0.005.")
    parser.add_argument("--plateau_tol", type=float, default=5.0, metavar="PERCENT",
                        help="error-plateau detection: stop if an iteration reduced the "
                             "error by less than PERCENT %% relative to the previous one "
                             "(default: 5). Only active from the 4th refinement step "
                             "onward. Set to 0 to disable the check.")
    parser.add_argument("--max_iter", type=int, default=10,
                        help="Maximum refinement iterations. Default: 10.")
    parser.add_argument("--refinement_factor", type=int, default=3,
                        help="Subdivisions per axis for refined regions. Default: 3.")
    parser.add_argument("--energy_window", type=float, default=0.1,
                        help="Ceiling for the mu-centred hotspot window [eV]. "
                             "Default: 0.1.")
    parser.add_argument("--keep_intermediate", action="store_true",
                        help="Keep intermediate mesh and output files instead of "
                             "deleting them.")
    parser.add_argument("--workdir", type=Path, default=Path.cwd(),
                        help="Directory for intermediate and output files. Default: cwd.")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable debug logging.")

    # WannGenerator options
    parser.add_argument("--charge", type=float, default=None,
                        help="Number of electrons in the projected bands. If omitted, "
                             "the charge is determined at the supplied chemical "
                             "potential and rounded, as in a plain lwann run.")
    parser.add_argument("--soc", action="store_true",
                        help="Wannier90 calculation includes spin-orbit coupling "
                             "(needed only for spin-unpolarised calculations).")
    parser.add_argument("--nocorrection", action="store_true",
                        help="Do not calculate the correction terms to the Peierls "
                             "approximation.")
    parser.add_argument("--intraonly", action="store_true",
                        help="Discard inter-band elements before writing. Also "
                             "avoids allocating the full (nband x nband) moment "
                             "arrays, which dominate memory for many-orbital "
                             "models.")
    parser.add_argument("--kblock", type=int, default=None, metavar="N",
                        help="k-points evaluated per pass when building H(k). "
                             "Default: sized automatically from --max-memory. "
                             "Results do not depend on this value (up to the "
                             "eigenvector phase of off-diagonal band elements).")
    parser.add_argument("--max-memory", dest="max_memory", type=float,
                        default=None, metavar="GB",
                        help="Budget [GB] for the transient per-block arrays. "
                             "Default: a quarter of the detected available "
                             "memory. Does not bound the output arrays, whose "
                             "size is reported at the start of each pass.")
    add_openmp_argument(parser)

    return parser.parse_args(argv)


def validate_inputs(args: argparse.Namespace) -> None:
    """Validate file paths and numeric ranges before starting the loop."""
    if args.initial_hdf5.suffix.lower() != ".hdf5":
        raise ValueError("Initial file must be an HDF5 (.hdf5) file.")
    if not args.initial_hdf5.exists():
        raise FileNotFoundError(args.initial_hdf5)
    if not args.wannier_dir.exists():
        raise FileNotFoundError(args.wannier_dir)
    if not args.wannier_dir.is_dir():
        raise ValueError(
            "wannier_dir must be the folder of a Wannier90 calculation, not a "
            "single file: {}".format(args.wannier_dir)
        )

    if args.error_tol <= 0:
        raise ValueError("error_tol must be positive.")
    if args.max_iter <= 0:
        raise ValueError("max_iter must be positive.")
    if args.refinement_factor < 1:
        raise ValueError("refinement_factor must be at least 1.")
    if args.kblock is not None and args.kblock < 1:
        raise ValueError("kblock must be at least 1.")
    if args.max_memory is not None and args.max_memory <= 0:
        raise ValueError("max-memory must be positive.")
    if args.charge is not None and args.charge < 0:
        raise ValueError("charge must be non-negative.")
    if args.T_min < 0:
        raise ValueError("Temperature must be non-negative")
    if args.T_min == 0.0:
        logger.warning(
            "T_min = 0 is not supported (beta = 1/(kB*T) diverges). "
            "Increasing to T_min = %g K.", T_MIN_FLOOR_K
            )
        args.T_min = T_MIN_FLOOR_K


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format=LOG_FORMAT,
    )

    logger.info(openmp_info_line(_openmp_ncores))

    try:
        validate_inputs(args)
    except Exception as exc:
        logger.error("Input validation failed: %s", exc)
        return 1

    if args.charge is not None:
        logger.warning(
            "--charge given: the chemical potential is recomputed from the "
            "charge on every mesh and may drift from the value %.6f used for "
            "hotspot detection. Omit --charge to keep mu fixed instead.",
            args.chemical_potential,
        )

    # Inherit the dimensionality of the coarse calculation rather than
    # re-deriving it from the refined k-points: a layered system can have a
    # 3D k-mesh while the Wannier projection is effectively 2D.
    source_dims = read_source_dims(args.initial_hdf5)

    # An irreducible source mesh covers only a wedge; the moments must then be
    # group-averaged over the star of every k-point (see _setCustomSymmetries).
    source_irreducible, source_symop = read_source_symmetry(args.initial_hdf5)
    if source_irreducible and source_symop is not None:
        logger.info("Source mesh is irreducible; refined meshes will be "
                    "symmetrised over %d operations.", source_symop.shape[0])

    generator = WannGenerator(
        directory         = args.wannier_dir,
        charge            = args.charge,
        mu                = args.chemical_potential,
        soc               = args.soc,
        peierlscorrection = not args.nocorrection,
        intraonly         = args.intraonly,
        dims              = source_dims,
        kblock            = args.kblock,
        memory_budget_gb  = args.max_memory,
        symop             = source_symop,
    )

    params = RefinementParams(
        initial_hdf5       = args.initial_hdf5,
        chemical_potential = args.chemical_potential,
        gamma_min          = args.gamma_min,
        T_min              = args.T_min,
        error_tol          = args.error_tol,
        max_iter           = args.max_iter,
        refinement_factor  = args.refinement_factor,
        energy_window      = args.energy_window,
        workdir            = args.workdir,
        keep_intermediate  = args.keep_intermediate,
        plateau_tol        = args.plateau_tol / 100.0,  # CLI takes percent
    )

    return run_refinement(params, generator)


if __name__ == "__main__":
    sys.exit(main())
