#!/usr/bin/env python3
"""
scripts/ltb/supercell.py  --  Tight-Binding Supercell Generator

Create a periodic or finite supercell from a tight-binding input file.

Called via:
  ltb supercell in.tbdata out.tbdata N1 N2 N3 [options]   # through the wrapper
  python -m scripts.ltb.supercell in.tbdata out.tbdata 2 2 1 [options]

Based on ltb_supercell (v0.1, J.M. Tomczak, 2026).
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
from collections import defaultdict

import numpy as np


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def read_tb_file(filename: str):
    lattice = []; atoms = []; orbitals = {}; hoppings = []; section = None

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.lstrip().startswith("#"):
                continue
            if line.startswith("begin"):
                section = line.split()[1]; continue
            if line.startswith("end"):
                section = None; continue
            line = line.split("#")[0].strip()
            if not line:
                continue
            parts = line.split()

            if section == "real_lattice":
                lattice.append(list(map(float, parts)))
            elif section == "atoms":
                atoms.append((int(parts[0]),
                              float(parts[1]), float(parts[2]), float(parts[3])))
            elif section == "orbitals":
                idx = int(parts[0])
                orbitals[idx] = np.array(list(map(float, parts[1:4])))
            elif section == "hopping":
                dx, dy, dz = map(int, parts[0:3])
                o1, o2     = map(int, parts[3:5])
                tr = float(parts[5])
                ti = float(parts[6]) if len(parts) == 7 else 0.0
                hoppings.append((dx, dy, dz, o1, o2, tr, ti))

    return np.array(lattice), atoms, orbitals, hoppings


def write_tb_file(filename: str, lattice, atoms, orbitals, hoppings) -> None:
    with open(filename, "w") as f:
        f.write("begin hopping\n")
        for dx, dy, dz, o1, o2, tr, ti in hoppings:
            if ti != 0.0:
                f.write(f"{dx} {dy} {dz} {o1} {o2} {tr} {ti}\n")
            else:
                f.write(f"{dx} {dy} {dz} {o1} {o2} {tr}\n")
        f.write("end hopping\n\n")

        f.write("begin atoms\n")
        for a in atoms:
            f.write(f"{a[0]} {a[1]} {a[2]} {a[3]}\n")
        f.write("end atoms\n\n")

        f.write("begin orbitals\n")
        for o in orbitals:
            f.write(f"{o[0]} {o[1]} {o[2]} {o[3]}\n")
        f.write("end orbitals\n\n")

        f.write("begin real_lattice\n")
        for v in lattice:
            f.write(f"{v[0]} {v[1]} {v[2]}\n")
        f.write("end real_lattice\n")


# ---------------------------------------------------------------------------
# Hermiticity check
# ---------------------------------------------------------------------------

def check_hermiticity(hoppings, label: str = "") -> int:
    """
    Verify that h(R, i, j) = h(-R, j, i)* for all entries.
    Prints a summary and returns the number of violations found.
    """
    hop_dict: dict = defaultdict(list)
    for dx, dy, dz, i, j, tr, ti in hoppings:
        hop_dict[(dx, dy, dz, i, j)].append(complex(tr, ti))

    errors = 0
    for dx, dy, dz, i, j, tr, ti in hoppings:
        key_conj = (-dx, -dy, -dz, j, i)
        val = complex(tr, ti)
        if key_conj not in hop_dict:
            print(f"  Missing conjugate partner for ({dx},{dy},{dz},{i},{j})")
            errors += 1
        else:
            if not any(np.isclose(val, np.conj(v)) for v in hop_dict[key_conj]):
                print(f"  Conjugate mismatch for ({dx},{dy},{dz},{i},{j}): "
                      f"{val} vs conj({hop_dict[key_conj]})")
                errors += 1

    tag = f" ({label})" if label else ""
    if errors == 0:
        print(f"Hermiticity check{tag}: OK")
    else:
        print(f"Hermiticity check{tag}: {errors} issue(s) found")
    return errors


# ---------------------------------------------------------------------------
# Supercell construction
# ---------------------------------------------------------------------------

def build_supercell(lattice, atoms, orbitals, hoppings,
                    N1: int, N2: int, N3: int,
                    make_finite: bool = False):
    """
    Expand a primitive tight-binding model into an N1 x N2 x N3 supercell.

    When no explicit 'begin orbitals' section is present in the TB file,
    the atoms list is used as a 1-to-1 fallback (one orbital per atom,
    indexed 1..N_atoms).

    Parameters
    ----------
    lattice      : (3,3) array of primitive lattice vectors.
    atoms        : list of (sort, x, y, z) in fractional coordinates.
    orbitals     : dict {idx -> (x,y,z)} in fractional coordinates.
                   May be empty if the TB file has no 'begin orbitals' block.
    hoppings     : list of (dx, dy, dz, o1, o2, t_real, t_imag).
    N1,N2,N3     : supercell repetitions along a1, a2, a3.
    make_finite  : if True, discard hoppings that cross the supercell boundary.

    Returns
    -------
    new_lattice, new_atoms, new_orbitals, new_hoppings
    """
    # Fall back to atoms when no explicit orbitals section is present.
    if not orbitals:
        orbitals = {i + 1: np.array([ax, ay, az])
                    for i, (_, ax, ay, az) in enumerate(atoms)}

    index_map: dict = {}
    new_orbitals = []
    new_atoms    = []
    counter      = 1

    for n1 in range(N1):
        for n2 in range(N2):
            for n3 in range(N3):
                for (sort_label, ax, ay, az), (old_idx, opos) in \
                        zip(atoms, orbitals.items()):
                    shift   = np.array([n1, n2, n3])
                    new_pos = (opos + shift) / np.array([N1, N2, N3])
                    new_orbitals.append((counter, *new_pos))
                    new_atoms.append((sort_label, *new_pos))
                    index_map[(n1, n2, n3, old_idx)] = counter
                    counter += 1

    new_hoppings = []
    for n1 in range(N1):
        for n2 in range(N2):
            for n3 in range(N3):
                for dx, dy, dz, o1, o2, tr, ti in hoppings:
                    tx = n1 + dx;  ty = n2 + dy;  tz = n3 + dz
                    SX = tx // N1; SY = ty // N2; SZ = tz // N3
                    wx = tx %  N1; wy = ty %  N2; wz = tz %  N3

                    if make_finite and (SX != 0 or SY != 0 or SZ != 0):
                        continue

                    i = index_map[(n1, n2, n3, o1)]
                    j = index_map[(wx, wy, wz, o2)]
                    new_hoppings.append((SX, SY, SZ, i, j, tr, ti))

    new_lattice = np.array([
        lattice[0] * N1,
        lattice[1] * N2,
        lattice[2] * N3,
    ])
    return new_lattice, new_atoms, new_orbitals, new_hoppings


def remove_duplicates(hoppings):
    return list(set(hoppings))


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def get_parser():
    return argparse.ArgumentParser(
        prog='ltb supercell',
        description='Create a periodic or finite supercell from a tight-binding input file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Functionality
-------------
  Expands orbitals into an N1 x N2 x N3 supercell.
  Reconstructs the hopping list with updated orbital indices.
  Primitive hoppings crossing the supercell boundary are:
    - wrapped and stored with supercell displacement (periodic, default)
    - discarded with --make-finite (finite flake)

  Orbital indices in the output always start at 1 and increase sequentially.

Remarks
-------
  Not yet tested for multi-orbital atoms; the atom list may be incorrect
  in that case.

Usage via wrapper (recommended)
--------------------------------
  ltb supercell in.tbdata out.tbdata N1 N2 N3 [options]

Usage via direct call (expert)
-------------------------------
  python -m scripts.ltb.supercell in.tbdata out.tbdata N1 N2 N3 [options]
  python /path/to/scripts/ltb/supercell.py in.tbdata out.tbdata N1 N2 N3 [options]

Examples
--------
  ltb supercell graphene.tbdata graphene_2x2.tbdata 2 2 1
  ltb supercell graphene.tbdata flake.tbdata 4 4 1 --make-finite --check-hermiticity-out
""",
    )


def _add_arguments(parser):
    parser.add_argument("infile",  help="Input TB file.")
    parser.add_argument("outfile", help="Output TB file.")
    parser.add_argument("N1", type=int, help="Supercell size along lattice vector a1.")
    parser.add_argument("N2", type=int, help="Supercell size along lattice vector a2.")
    parser.add_argument("N3", type=int, help="Supercell size along lattice vector a3.")
    parser.add_argument("--make-finite", action="store_true",
                        help="Discard hoppings that cross the supercell boundary "
                             "(produces a finite flake with no periodic boundary conditions).")
    parser.add_argument("--remove-duplicates", action="store_true",
                        help="Remove identical hopping entries after construction.")
    parser.add_argument("--check-hermiticity-in", action="store_true",
                        help="Check hermiticity of the input hopping list.")
    parser.add_argument("--check-hermiticity-out", action="store_true",
                        help="Check hermiticity of the generated hopping list.")
    return parser


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = get_parser()
    _add_arguments(parser)
    args = parser.parse_args(argv)

    lattice, atoms, orbitals, hoppings = read_tb_file(args.infile)

    if args.check_hermiticity_in:
        check_hermiticity(hoppings, "input")

    new_lattice, new_atoms, new_orbitals, new_hoppings = build_supercell(
        lattice, atoms, orbitals, hoppings,
        args.N1, args.N2, args.N3,
        make_finite=args.make_finite,
    )

    if args.remove_duplicates:
        new_hoppings = remove_duplicates(new_hoppings)

    if args.check_hermiticity_out:
        check_hermiticity(new_hoppings, "output")

    write_tb_file(args.outfile, new_lattice, new_atoms, new_orbitals, new_hoppings)
    print(f"Supercell written to {args.outfile}  "
          f"({len(new_orbitals)} orbitals, {len(new_hoppings)} hoppings)")


if __name__ == "__main__":
    main()
