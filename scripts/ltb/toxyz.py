#!/usr/bin/env python3
"""
scripts/ltb/toxyz.py  --  Convert a LinReTraCe tight-binding file to XYZ format.

Features
--------
- Reads atomic positions from the 'atoms' section (fractional coordinates).
- Converts to Cartesian coordinates using the real_lattice vectors.
- Maps atom sort indices to chemical element symbols via --atom flags.
- Warns if periodic hoppings are present (they are irrelevant for XYZ output).
- Optional centering of the structure at the origin.
- Writes to stdout or a file.

Called via:
  ltb toxyz input.tbdata --atom 1 C [options]        # through the wrapper
  python -m scripts.ltb.toxyz input.tbdata --atom 1 C
  python /path/to/scripts/ltb/toxyz.py input.tbdata --atom 1 C

Based on tb2xyz (J.M. Tomczak, 2026).

Notes on atoms vs orbitals
---------------------------
The 'atoms' section defines atomic positions used for XYZ output and
visualisation.  The 'orbitals' section defines the positions at which
orbitals sit for the generalised Peierls substitution (vector potential
coupling).  For XYZ output we always use 'atoms', since XYZ is an
atomic-position format.  The 'orbitals' section is ignored here.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure package root is on sys.path when invoked directly by path.
# scripts/ltb/ is two levels below the root.
_root = Path(__file__).resolve().parents[2]
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

import argparse

import numpy as np


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_tbdata(filename: str):
    """
    Read atoms, lattice vectors, and check for periodic hoppings.

    Returns
    -------
    atoms : list of (sort_index, frac_coords)
    lattice : (3,3) ndarray of row vectors [Angstrom]
    hopping_has_periodic : bool  -- True if any hopping has (Rx,Ry,Rz) != (0,0,0)
    """
    atoms: list = []
    lattice: list = []
    hopping_has_periodic = False
    section = None

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("begin"):
                section = line.split()[1]
                continue
            elif line.startswith("end"):
                section = None
                continue

            if section == "atoms":
                parts = line.split()
                atom_type = int(parts[0])
                frac = np.array(list(map(float, parts[1:4])))
                atoms.append((atom_type, frac))

            elif section == "real_lattice":
                lattice.append(list(map(float, line.split())))

            elif section == "hopping":
                parts = line.split()
                Rx, Ry, Rz = map(int, parts[0:3])
                if Rx != 0 or Ry != 0 or Rz != 0:
                    hopping_has_periodic = True

    if len(lattice) != 3:
        raise ValueError(
            f"real_lattice section must contain exactly 3 vectors, found {len(lattice)}."
        )

    return atoms, np.array(lattice), hopping_has_periodic


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def fractional_to_cartesian(frac: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """Convert fractional coordinates to Cartesian: cart = frac @ lattice."""
    return frac @ lattice


def center_structure(cart_atoms: list) -> list:
    """Shift all atoms so their centroid is at the origin."""
    coords = np.array([c for _, c in cart_atoms])
    centre = coords.mean(axis=0)
    return [(el, coord - centre) for el, coord in cart_atoms]


# ---------------------------------------------------------------------------
# XYZ output
# ---------------------------------------------------------------------------

def write_xyz(cart_atoms: list, comment: str, outfile: str | None = None) -> None:
    lines = [str(len(cart_atoms)), comment]
    for element, coord in cart_atoms:
        lines.append(f"{element} {coord[0]:.10f} {coord[1]:.10f} {coord[2]:.10f}")
    content = "\n".join(lines) + "\n"

    if outfile:
        with open(outfile, "w") as f:
            f.write(content)
    else:
        print(content, end="")


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def get_parser():
    return argparse.ArgumentParser(
        prog="ltb toxyz",
        description=(
            "Convert a LinReTraCe tight-binding file to XYZ format.\n\n"
            "Atom sort indices (from the 'atoms' section) are mapped to chemical\n"
            "element symbols via --atom flags.  Unmapped types are written as XN."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage via wrapper (recommended)
--------------------------------
  ltb toxyz input.tbdata --atom 1 C [options]
  ltb toxyz input.tbdata --atom 1 C --atom 2 H --center -o output.xyz

Usage via direct call (expert)
-------------------------------
  python -m scripts.ltb.toxyz input.tbdata --atom 1 C
  python /path/to/scripts/ltb/toxyz.py input.tbdata --atom 1 C

Examples
--------
  ltb toxyz graphene.tbdata --atom 1 C
  ltb toxyz BN.tbdata --atom 1 B --atom 2 N --center -o BN.xyz
""",
    )


def _add_arguments(parser):
    parser.add_argument("inputfile",
                        help="Input tight-binding file (.tbdata or similar).")
    parser.add_argument("-o", "--output", default=None, metavar="FILE",
                        help="Output XYZ file. Default: stdout.")
    parser.add_argument("--center", action="store_true",
                        help="Shift structure so its centroid is at the origin.")
    parser.add_argument("--atom", nargs=2, metavar=("TYPE", "ELEMENT"),
                        action="append", default=[],
                        help="Map atom sort index TYPE to chemical element ELEMENT. "
                             "Repeat for each type: --atom 1 C --atom 2 H. "
                             "Unmapped types are written as XN.")
    return parser


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = get_parser()
    _add_arguments(parser)
    args = parser.parse_args(argv)

    atom_map = {}
    for type_str, element in args.atom:
        try:
            atom_map[int(type_str)] = element
        except ValueError:
            parser.error(f"Atom type must be an integer, got '{type_str}'.")

    atoms, lattice, periodic = parse_tbdata(args.inputfile)

    if not atoms:
        sys.exit(
            "ltb toxyz: no atoms found in the input file.\n"
            "Make sure the file contains a 'begin atoms ... end atoms' section."
        )

    if periodic:
        print(
            "WARNING: non-zero lattice translations detected in the hopping section.\n"
            "These are irrelevant for XYZ output (only atomic positions are written).",
            file=sys.stderr,
        )

    unmapped = sorted({t for t, _ in atoms if t not in atom_map})
    if unmapped:
        print(
            f"WARNING: atom sort indices {unmapped} have no --atom mapping "
            f"and will be written as X1, X2, ... etc.",
            file=sys.stderr,
        )

    cart_atoms = []
    for atom_type, frac in atoms:
        cart = fractional_to_cartesian(frac, lattice)
        element = atom_map.get(atom_type, f"X{atom_type}")
        cart_atoms.append((element, cart))

    if args.center:
        cart_atoms = center_structure(cart_atoms)

    write_xyz(
        cart_atoms,
        comment=f"Converted from {args.inputfile}",
        outfile=args.output,
    )

    if args.output:
        print(f"XYZ written to {args.output}  ({len(cart_atoms)} atoms)")


if __name__ == "__main__":
    main()
