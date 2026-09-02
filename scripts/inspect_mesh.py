#!/usr/bin/env python3
"""
inspect-mesh: report what a k-mesh can and cannot resolve at a given kernel width.

Answers two questions that are NOT the same question and need two different
criteria:

  density      Is the mesh fine enough?  Judged by the sum-rule metric, which
               on a regular grid tracks the transport error monotonically.
  aspect ratio Are the axes divided in the right proportion?  The metric is
               BLIND to this -- it ranks a 32x32 mesh above a 64x8 one that is
               250x more accurate -- because mesh anisotropy is a k-space
               property while the metric lives on the energy axis.  Judged
               instead by the band-diagonal velocity tensor.

Usage
-----
    ltb   inspect-mesh mesh.hdf5 [--T K] [--gamma EV] [--target TOL]
    lwann inspect-mesh mesh.hdf5 [--T K] [--gamma EV] [--target TOL]

Give the WIDEST corner of the planned sweep: a mesh adequate at 1 meV can be
useless at 26 meV, and the wide corner is the one a starting mesh has to get
right, because refinement cannot repair it afterwards.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

logger = logging.getLogger("inspect-mesh")


def build_parser(prog: str = "inspect-mesh") -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Report the resolving power and axis proportions of a k-mesh.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("mesh", type=Path,
                        help="Energy/k-mesh HDF5 file to inspect.")
    parser.add_argument("--T", type=float, default=300.0, metavar="K",
                        help="Temperature of the corner to test [K]. Give the "
                             "WIDEST temperature of the planned sweep. Default: 300.")
    parser.add_argument("--gamma", type=float, default=0.01, metavar="EV",
                        help="Scattering rate of the corner to test [eV]. Give the "
                             "LARGEST gamma of the planned sweep; when gamma is "
                             "temperature dependent, its maximum over the sweep is "
                             "the safe choice, since that can only over-estimate the "
                             "kernel width, never under-estimate it. Default: 0.01.")
    parser.add_argument("--target", type=float, default=5e-3, metavar="TOL",
                        help="Metric target to judge the density against; the "
                             "divisions needed to reach it are then suggested. "
                             "Defaults to 5e-3, the same default as refine's "
                             "--error_tol and --parent_tol, so that inspect-mesh and "
                             "refine agree on what 'converged' means. Pass 0 to "
                             "report the metric without a verdict.")
    parser.add_argument("--debug", action="store_true", help="Verbose logging.")
    return parser


def main(argv=None) -> int:
    args = build_parser(sys.argv[0]).parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    # Imported here so that --help works without numpy/h5py present.
    from scripts.kmesh_refinement import format_inspection, inspect_mesh

    if not args.mesh.is_file():
        logger.error("No such file: %s", args.mesh)
        return 1

    target = None if args.target is not None and args.target <= 0 else args.target
    # NOTE: no --energy_window here.  That parameter enters hotspot SELECTION
    # during refinement, never compute_error/panel_defects, so it would have no
    # effect on anything this command reports.
    info = inspect_mesh(args.mesh, args.T, args.gamma, target=target)
    print(format_inspection(info))

    if target is not None and info["error"] > target:
        return 2          # distinguishable exit code for scripting
    return 0


if __name__ == "__main__":
    sys.exit(main())
