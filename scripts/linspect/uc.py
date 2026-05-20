#!/usr/bin/env python3
"""
scripts/linspect/uc.py  -  Visualise the real-space unit cell of a tight-binding system.

Currently supports 2D systems only (z-component of hoppings is ignored).
3D systems are detected automatically and rejected with a clear error message.

Called via:
  linspect uc graphene.tbdata [options]           # through the wrapper
  python -m scripts.linspect.uc graphene.tbdata   # direct module call
  python /path/to/scripts/linspect/uc.py ...      # direct path call (expert)

Based on ltb_show2d (v0.3, J.M. Tomczak, 2026).
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure package root is on sys.path when invoked directly by path.
# scripts/linspect/ is two levels below the root.
_root = Path(__file__).resolve().parents[2]
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

import argparse
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def read_tb_file(filename):
    lattice_vectors = np.eye(3)
    atoms           = {}
    orbitals        = {}
    hoppings        = []
    onsite          = defaultdict(float)
    section         = None
    lattice_lines   = []

    with open(filename) as f:
        for line in f:
            line = line.split("#")[0].strip()
            if not line:
                continue
            if line.startswith("begin"):
                section = line.split()[1]
                continue
            elif line.startswith("end"):
                if section == "real_lattice":
                    lattice_vectors = np.array(lattice_lines)
                section = None
                lattice_lines = []
                continue

            parts = line.split()
            if section == "real_lattice":
                lattice_lines.append(list(map(float, parts)))
            elif section == "atoms":
                idx = len(atoms) + 1
                pos = np.array(list(map(float, parts[1:4])))
                atoms[idx] = pos
            elif section == "orbitals":
                idx = len(orbitals) + 1
                pos = np.array(list(map(float, parts[1:4])))
                orbitals[idx] = pos
            elif section == "hopping":
                dx, dy, dz = map(int, parts[0:3])
                o1, o2     = int(parts[3]), int(parts[4])
                hval       = float(parts[5])
                if o1 == o2 and dx == dy == dz == 0:
                    onsite[o1] = hval
                else:
                    hoppings.append((dx, dy, dz, o1, o2, hval))

    return lattice_vectors, atoms, orbitals, hoppings, onsite


# ---------------------------------------------------------------------------
# Dimension check
# ---------------------------------------------------------------------------

def check_2d(hoppings, lattice_vectors):
    """
    Abort with a helpful message if the system appears to be 3D.

    A system is treated as 3D if any hopping has dz != 0, or if the
    z-lattice vector has a significant out-of-plane component.
    """
    if any(h[2] != 0 for h in hoppings):
        sys.exit(
            "linspect uc: this system has hoppings along z (dz != 0) and "
            "appears to be 3D.\n"
            "3D unit-cell visualisation is not yet supported.\n"
            "If your system is truly 2D, check that dz=0 for all hopping lines."
        )

    lz = np.linalg.norm(lattice_vectors[2, :2])   # in-plane component of c
    cz = abs(lattice_vectors[2, 2])                # out-of-plane component of c
    if cz < 1e-6 and lz > 1e-6:
        sys.exit(
            "linspect uc: the third lattice vector has no significant z-component "
            "and appears to describe a 3D direction in-plane.\n"
            "3D unit-cell visualisation is not yet supported."
        )


# ---------------------------------------------------------------------------
# Topology detection
# ---------------------------------------------------------------------------

def detect_topology(hoppings):
    """
    Returns one of:
        'finite'      - all hoppings within cell (0,0,*)
        'semi_x'      - periodic only along x  (dx!=0, all dy==0)
        'semi_y'      - periodic only along y  (dy!=0, all dx==0)
        'periodic'    - periodic in both x and y
    """
    has_dx = any(h[0] != 0 for h in hoppings)
    has_dy = any(h[1] != 0 for h in hoppings)

    if not has_dx and not has_dy:
        return 'finite'
    elif has_dx and not has_dy:
        return 'semi_x'
    elif has_dy and not has_dx:
        return 'semi_y'
    else:
        return 'periodic'


# ---------------------------------------------------------------------------
# Build real-space patch
# ---------------------------------------------------------------------------

def make_finite_patch(lattice, positions_dict, hoppings, Nx=4, Ny=4):
    """
    Build a finite real-space patch appropriate for the detected topology.

    Padding cells are added beyond the display region so that every hopping
    originating inside the display region has a valid target, giving every
    displayed atom a complete, uniform set of hopping lines.  Padding cells
    appear as hopping endpoints but are NOT rendered as atoms.

    Returns
    -------
    pos_dict            : {global_idx -> 3-vector (A)}  display + padding sites
    draw_sites          : set of global_idx to render as atoms (display region only)
    hoppings_finite     : [(i, j, hval), ...]
    topology            : 'finite' | 'semi_x' | 'semi_y' | 'periodic'
    primitive_of_global : {global_idx -> orbital_idx}
    """
    topology = detect_topology(hoppings)

    if topology == 'finite':
        nx_tile, ny_tile = 1, 1
        periodic_x, periodic_y = False, False
    elif topology == 'semi_x':
        nx_tile, ny_tile = Nx, 1
        periodic_x, periodic_y = True, False
    elif topology == 'semi_y':
        nx_tile, ny_tile = 1, Ny
        periodic_x, periodic_y = False, True
    else:  # 'periodic'
        nx_tile, ny_tile = Nx, Ny
        periodic_x, periodic_y = True, True

    pad_x_lo = max(0, -min((h[0] for h in hoppings), default=0)) if periodic_x else 0
    pad_x_hi = max(0,  max((h[0] for h in hoppings), default=0)) if periodic_x else 0
    pad_y_lo = max(0, -min((h[1] for h in hoppings), default=0)) if periodic_y else 0
    pad_y_hi = max(0,  max((h[1] for h in hoppings), default=0)) if periodic_y else 0

    nx_eff = nx_tile + pad_x_lo + pad_x_hi
    ny_eff = ny_tile + pad_y_lo + pad_y_hi

    primitive_of_global = {}
    pos_dict   = {}
    index_map  = {}
    global_idx = 1

    for ix in range(nx_eff):
        for iy in range(ny_eff):
            R = ix * lattice[0] + iy * lattice[1]
            for orb_idx, pos_frac in positions_dict.items():
                pos_real                   = pos_frac @ lattice + R
                pos_dict[global_idx]       = pos_real
                index_map[(ix, iy, orb_idx)] = global_idx
                primitive_of_global[global_idx] = orb_idx
                global_idx += 1

    draw_sites = set()
    for ix in range(pad_x_lo, pad_x_lo + nx_tile):
        for iy in range(pad_y_lo, pad_y_lo + ny_tile):
            for orb_idx in positions_dict:
                draw_sites.add(index_map[(ix, iy, orb_idx)])

    hoppings_finite = []
    for ix in range(pad_x_lo, pad_x_lo + nx_tile):
        for iy in range(pad_y_lo, pad_y_lo + ny_tile):
            for dx, dy, dz, o1, o2, hval in hoppings:
                key_from = (ix,      iy,      o1)
                key_to   = (ix + dx, iy + dy, o2)
                if key_from in index_map and key_to in index_map:
                    gi = index_map[key_from]
                    gj = index_map[key_to]
                    if gi in draw_sites and gj in draw_sites:
                        hoppings_finite.append((gi, gj, hval))

    R_shift  = pad_x_lo * lattice[0] + pad_y_lo * lattice[1]
    pos_dict = {k: v - R_shift for k, v in pos_dict.items() if k in draw_sites}

    return pos_dict, draw_sites, hoppings_finite, topology, primitive_of_global


# ---------------------------------------------------------------------------
# Unit-cell drawing
# ---------------------------------------------------------------------------

def _draw_unitcell(ax, lattice, topology, pos_dict, draw_sites, primitive_of_global):
    """
    Draw a red dashed parallelogram indicating the primitive unit cell.

    For 'periodic': the cell is drawn at the origin corner of the patch.
    For 'semi_x' / 'semi_y': the cell spans the full periodic direction and
    is centred on the actual atom extent in the finite direction.
    """
    a1 = lattice[0, :2]
    a2 = lattice[1, :2]

    seen      = set()
    prim_cart = []
    for gi in sorted(draw_sites):
        orb = primitive_of_global[gi]
        if orb not in seen:
            prim_cart.append(pos_dict[gi][:2])
            seen.add(orb)
    centroid = np.array(prim_cart).mean(axis=0)

    cell_centre = 0.5 * (a1 + a2)
    shift       = centroid - cell_centre

    mat = np.column_stack([a1, a2])
    if abs(np.linalg.det(mat)) > 1e-10:
        frac  = np.linalg.solve(mat, shift)
        frac  = frac - np.round(frac)
        shift = mat @ frac

    if topology in ('periodic', 'semi_x', 'semi_y'):
        corners = np.array([[0, 0], a1, a1 + a2, a2, [0, 0]]) + shift
        ax.plot(corners[:, 0], corners[:, 1],
                linestyle='--', color='red', linewidth=2, zorder=4,
                label='unit cell')


# ---------------------------------------------------------------------------
# Main plotting routine
# ---------------------------------------------------------------------------

def plot_uc(filename, Nx=4, Ny=4,
            show_hoppings=True, show_onsite=True, show_unitcell=True,
            cmap='viridis', xlim=None, ylim=None, out=None):

    lattice, atoms, orbitals, hoppings, onsite = read_tb_file(filename)
    check_2d(hoppings, lattice)
    positions_dict = orbitals if orbitals else atoms

    pos_dict, draw_sites, hoppings_finite, topology, primitive_of_global = \
        make_finite_patch(lattice, positions_dict, hoppings, Nx=Nx, Ny=Ny)

    print(f"[linspect uc] Detected topology: {topology}")

    sorted_draw  = sorted(draw_sites)
    positions_2d = np.array([pos_dict[idx][:2] for idx in sorted_draw])
    all_pos_2d   = {idx: pos_dict[idx][:2] for idx in pos_dict}

    if show_onsite:
        energies = np.array([
            onsite.get(primitive_of_global[idx], 0.0)
            for idx in sorted_draw
        ])
        norm   = plt.Normalize(vmin=np.min(energies), vmax=np.max(energies))
        colors = plt.colormaps[cmap](norm(energies))
    else:
        colors = 'orange'

    fig, ax = plt.subplots(figsize=(8, 6))

    if show_hoppings and hoppings_finite:
        hvals        = np.array([abs(h[2]) for h in hoppings_finite])
        hmin, hmax   = hvals.min(), hvals.max()
        for i, j, hval in hoppings_finite:
            lw = 1 + 2 * (abs(hval) - hmin) / (hmax - hmin) if hmax > hmin else 1
            p1 = all_pos_2d[i]
            p2 = all_pos_2d[j]
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 'gray', lw=lw, zorder=1)

    ax.scatter(positions_2d[:, 0], positions_2d[:, 1],
               c=colors, s=80, edgecolor='black', linewidth=1.2, zorder=3)

    if topology != 'finite' and show_unitcell:
        _draw_unitcell(ax, lattice, topology, pos_dict, draw_sites, primitive_of_global)

    pad = 0.1
    if xlim is not None:
        ax.set_xlim(xlim[0] - pad, xlim[1] + pad)
    if ylim is not None:
        ax.set_ylim(ylim[0] - pad, ylim[1] + pad)

    ax.set_aspect('equal')
    ax.set_xlabel('x [Å]')
    ax.set_ylabel('y [Å]')
    topology_label = {
        'finite':   'Finite cluster',
        'semi_x':   'Semi-infinite (periodic x, finite y)',
        'semi_y':   'Semi-infinite (periodic y, finite x)',
        'periodic': 'Periodic bulk',
    }[topology]
    ax.set_title(f'Tight-binding unit cell (2D) — {topology_label}')

    if show_onsite:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label('Onsite potential [eV]')

    plt.tight_layout()

    if out is not None:
        if '.' not in Path(out).name:
            out += '.pdf'
        fmt = out.split('.')[-1].lower()
        if fmt not in ('png', 'pdf', 'svg'):
            print(f"Warning: unknown extension '{fmt}', defaulting to png.")
            fmt = 'png'
        plt.savefig(out, dpi=300, format=fmt, bbox_inches='tight')
        print(f'Figure saved to {out}')
    else:
        plt.show()


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def get_parser():
    parser = argparse.ArgumentParser(
        prog='linspect uc',
        description='Visualise the real-space unit cell of a 2D tight-binding system.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
  linspect uc graphene.tbdata
  linspect uc graphene.tbdata --Nx 4 --Ny 4
  linspect uc square_flake.tbdata --no-hoppings
  linspect uc cuprate.tbdata --Nx 4 --Ny 4 --xlim 2 5 --ylim 2 5 --no-unitcell
  linspect uc slab.tbdata --Ny 6 --out slab.pdf
""")
    parser.add_argument('filename',
                        help='Tight-binding input file (.tbdata or similar).')
    parser.add_argument('--Nx', type=int, default=4,
                        help='Unit cells to replicate along x (periodic/semi-x). Default: 4.')
    parser.add_argument('--Ny', type=int, default=4,
                        help='Unit cells to replicate along y (periodic/semi-y). Default: 4.')
    parser.add_argument('--no-hoppings', action='store_true',
                        help='Disable drawing of hopping lines.')
    parser.add_argument('--no-onsite', action='store_true',
                        help='Disable onsite-energy colour coding.')
    parser.add_argument('--no-unitcell', action='store_true',
                        help='Disable unit-cell parallelogram (periodic/semi-infinite).')
    parser.add_argument('--cmap', default='viridis',
                        help='Matplotlib colormap for onsite energies. Default: viridis.')
    parser.add_argument('--xlim', nargs=2, type=float, metavar=('XMIN', 'XMAX'),
                        help='Manual x-axis limits [Å].')
    parser.add_argument('--ylim', nargs=2, type=float, metavar=('YMIN', 'YMAX'),
                        help='Manual y-axis limits [Å].')
    parser.add_argument('--out', metavar='FILE',
                        help='Save figure to FILE (png/pdf/svg) instead of showing.')
    return parser


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main(argv=None):
    args = get_parser().parse_args(argv)
    plot_uc(
        args.filename,
        Nx             = args.Nx,
        Ny             = args.Ny,
        show_hoppings  = not args.no_hoppings,
        show_onsite    = not args.no_onsite,
        show_unitcell  = not args.no_unitcell,
        cmap           = args.cmap,
        xlim           = args.xlim,
        ylim           = args.ylim,
        out            = args.out,
    )


if __name__ == '__main__':
    main()
