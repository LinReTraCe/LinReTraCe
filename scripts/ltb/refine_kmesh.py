#!/usr/bin/env python3
"""
scripts/ltb/refine_kmesh.py

ltb entry point for iterative k-mesh refinement.

Parses command-line arguments, constructs an LtbGenerator, translates the
arguments into a RefinementParams object, and hands both to the shared
generator-agnostic loop in scripts/refine_kmesh.py.

Usage
-----
    python -m scripts.ltb.refine_kmesh initial_hdf5 tb_file chemical_potential
                                        gamma_min t_min [options]

Positional arguments
--------------------
    initial_hdf5        Initial tight-binding HDF5 file (must end with .hdf5).
    tb_file             Wannier-style tight-binding input file.
    chemical_potential  Chemical potential (mu) used for df/da and the generator.
    gamma_min           Minimum gamma parameter for df/da.
    t_min               Minimum temperature for df/da.

Optional arguments
------------------
    --error_tol ERR         Target df/da error threshold (default: 0.005).
    --max_iter N            Maximum refinement iterations (default: 10).
    --refinement_factor F   Subdivisions per axis for refined hotspots (default: 3).
    --energy_window E       mu-centred energy window for hotspot detection (default: 0.1).
    --filling FILL          Electron filling for the TB calculation (default: 2.0).
    --mushift               Shift energies so that mu = 0 after diagonalisation.
    --corronly              Use only the multi-orbital Peierls correction term.
    --vector                Also save full Hamiltonian matrices to the output HDF5.
    --intra VALUE           Set all intra-band optical elements to VALUE.
    --inter VALUE           Set all inter-band optical elements to VALUE.
    --intraonly             Discard inter-band elements before writing.
    --keep_intermediate     Keep intermediate mesh and HDF5 files.
    --workdir DIR           Directory for intermediate and output files (default: cwd).
    --verbose               Enable debug logging.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure the package root is on sys.path when this script is run directly
# by absolute path. scripts/ltb/ sits two levels below the root.
_root = Path(__file__).resolve().parents[2]
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

import argparse
import logging


from structure.generators.ltb_gen import LtbGenerator
from scripts.kmesh_refinement import RefinementParams, run_refinement

LOG_FORMAT = "%(levelname)s: %(message)s"
logger = logging.getLogger(__name__)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog='ltb refine',
        description='Iteratively refine a LinReTraCe k-mesh towards a target error threshold.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage via wrapper (recommended)
--------------------------------
  ltb refine initial.hdf5 tb_file mu gamma_min t_min [options]

Usage via direct call (expert)
--------------------------------
  python -m scripts.ltb.refine_kmesh initial.hdf5 tb_file mu gamma_min t_min [options]
  python /path/to/scripts/ltb/refine_kmesh.py initial.hdf5 tb_file mu gamma_min t_min [options]

Positional arguments
--------------------
  initial_hdf5        Initial tight-binding HDF5 file (must end with .hdf5).
  tb_file             Wannier-style tight-binding input file.
  chemical_potential  Chemical potential (mu) used for df/da and the generator.
  gamma_min           Minimum gamma parameter for df/da.
  t_min               Minimum temperature for df/da.
""",
    )

    # Positional arguments -- identical to the original ltb_refine_kmesh,
    # except --ltb_use_custom_kmesh is gone (no longer needed).
    parser.add_argument("initial_hdf5", type=Path,
                        help="Initial tight-binding HDF5 file (must end with .hdf5).")
    parser.add_argument("tb_file", type=Path,
                        help="Wannier-style tight-binding input file.")
    parser.add_argument("chemical_potential", type=float,
                        help="Chemical potential used for refinement and the TB calculation.")
    parser.add_argument("gamma_min", type=float,
                        help="Minimum gamma parameter for df/da.")
    parser.add_argument("t_min", type=float,
                        help="Minimum temperature for df/da.")

    # Refinement loop options
    parser.add_argument("--error_tol", type=float, default=5e-3,
                        help="Target error threshold. Default: 0.005.")
    parser.add_argument("--max_iter", type=int, default=10,
                        help="Maximum refinement iterations. Default: 10.")
    parser.add_argument("--refinement_factor", type=int, default=3,
                        help="Subdivisions per axis for refined regions. Default: 3.")
    parser.add_argument("--energy_window", type=float, default=0.1,
                        help="Energy window around mu for hotspot detection. Default: 0.1.")
    parser.add_argument("--keep_intermediate", action="store_true",
                        help="Keep intermediate mesh and output files instead of overwriting them.")
    parser.add_argument("--workdir", type=Path, default=Path.cwd(),
                        help="Directory for intermediate and output files. Default: cwd.")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable debug logging.")

    # LtbGenerator options
    parser.add_argument("--filling", type=float, default=2.0,
                        help="Electron filling for the TB calculation. Default: 2.0.")
    parser.add_argument("--mushift", action="store_true",
                        help="Shift energies so that mu = 0 after diagonalisation.")
    parser.add_argument("--corronly", action="store_true",
                        help="Use only the multi-orbital Peierls correction term.")
    parser.add_argument("--vector", action="store_true",
                        help="Also save full Hamiltonian matrices to the output HDF5.")
    parser.add_argument("--intra", type=float, default=None,
                        help="Set all intra-band optical elements to VALUE.")
    parser.add_argument("--inter", type=float, default=None,
                        help="Set all inter-band optical elements to VALUE.")
    parser.add_argument("--intraonly", action="store_true",
                        help="Discard inter-band elements before writing.")

    return parser.parse_args(argv)


def validate_inputs(args: argparse.Namespace) -> None:
    """Validate file paths and numeric ranges before starting the loop."""
    if args.initial_hdf5.suffix.lower() != ".hdf5":
        raise ValueError("Initial file must be an HDF5 (.hdf5) file.")
    if not args.initial_hdf5.exists():
        raise FileNotFoundError(args.initial_hdf5)
    if not args.tb_file.exists():
        raise FileNotFoundError(args.tb_file)

    if args.error_tol <= 0:
        raise ValueError("error_tol must be positive.")
    if args.max_iter <= 0:
        raise ValueError("max_iter must be positive.")
    if args.refinement_factor < 1:
        raise ValueError("refinement_factor must be at least 1.")
    if args.t_min <= 0:
        raise ValueError("Temperature must be positive for df/da calculation.")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format=LOG_FORMAT,
    )

    try:
        validate_inputs(args)
    except Exception as exc:
        logger.error("Input validation failed: %s", exc)
        return 1

    generator = LtbGenerator(
        tb_file   = args.tb_file,
        filling   = args.filling,
        mu        = args.chemical_potential,
        mushift   = args.mushift,
        corronly  = args.corronly,
        vector    = args.vector,
        intra     = args.intra,
        inter     = args.inter,
        intraonly = args.intraonly,
    )

    params = RefinementParams(
        initial_hdf5       = args.initial_hdf5,
        chemical_potential = args.chemical_potential,
        gamma_min          = args.gamma_min,
        t_min              = args.t_min,
        error_tol          = args.error_tol,
        max_iter           = args.max_iter,
        refinement_factor  = args.refinement_factor,
        energy_window      = args.energy_window,
        workdir            = args.workdir,
        keep_intermediate  = args.keep_intermediate,
    )

    return run_refinement(params, generator)


if __name__ == "__main__":
    sys.exit(main())
