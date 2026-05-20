#!/usr/bin/env python3
"""
scripts/linspect/bz.py  -  Visualise the Brillouin-zone k-mesh of a LinReTraCe HDF5 file.

Called via:
  linspect bz run.hdf5 [options]           # through the wrapper
  python -m scripts.linspect.bz run.hdf5  # direct module call
  python /path/to/scripts/linspect/bz.py  # direct path call (expert)
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

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.spatial import Voronoi


# ------------------------------------------------------------------------------
# BZ construction via Voronoi
# ------------------------------------------------------------------------------

def bz_polygon_2d(kvec_rows, ax1, ax2):
    """
    Return the first Brillouin zone as a closed (N,2) array of Cartesian
    vertices in the plane spanned by axes ax1 and ax2.

    The BZ is the Wigner-Seitz cell of the reciprocal lattice, computed as
    the Voronoi cell around the origin among a 9x9x9 supercell of G-vectors.
    """
    ns = range(-4, 5)
    G_list = []
    for n1 in ns:
        for n2 in ns:
            for n3 in ns:
                G_list.append(n1*kvec_rows[0] + n2*kvec_rows[1] + n3*kvec_rows[2])
    G    = np.array(G_list)
    G2d  = G[:, [ax1, ax2]]
    origin_idx = np.argmin(np.linalg.norm(G2d, axis=1))

    vor = Voronoi(G2d)
    region_idx    = vor.point_region[origin_idx]
    vertex_indices = vor.regions[region_idx]

    if -1 in vertex_indices:
        raise RuntimeError("BZ Voronoi cell is open; try a denser G grid.")

    verts  = vor.vertices[vertex_indices]
    centre = verts.mean(axis=0)
    angles = np.arctan2(verts[:, 1] - centre[1], verts[:, 0] - centre[0])
    verts  = verts[np.argsort(angles)]
    return np.vstack([verts, verts[0]])   # closed polygon


# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

AXIS_LABELS = {'x': r'$k_x$', 'y': r'$k_y$', 'z': r'$k_z$'}
AXIS_INDEX  = {'x': 0, 'y': 1, 'z': 2}


def frac_to_cart(kpoints_frac, kvec_rows):
    """Convert fractional k-coordinates to Cartesian: k_cart = k_frac @ kvec."""
    return kpoints_frac @ kvec_rows   # (N, 3)


def auto_scale(kpts_cart, window=27):
    """
    Estimate a marker scale factor from the k-point distribution.

    Strategy: sort k-points lexicographically in Cartesian space, then for
    each point find its nearest neighbour within a window of W adjacent
    entries in the sorted list.  This is O(N * W) — linear in N for fixed W.

    W=27 safely covers a 3x3x3 neighbourhood, which is enough to find true
    nearest neighbours even for non-orthogonal lattices and after multiple
    refinement iterations where coarse and fine points are interleaved.

    Simple consecutive differences (W=1) fail for non-orthogonal meshes
    because lexicographic sorting creates large row-crossing jumps that
    are not nearest-neighbour distances.

    Returns d_min / d_max where d_min and d_max are the global minimum and
    maximum nearest-neighbour distances across all k-points:
      - Uniform mesh:  d_min == d_max  ->  scale = 1.0  (no shrinkage)
      - Refined mesh:  d_min <  d_max  ->  scale < 1.0  (markers shrink)
    """
    N = kpts_cart.shape[0]
    if N < 2:
        return 1.0

    order     = np.lexsort(kpts_cart[:, ::-1].T)
    sorted_pts = kpts_cart[order]
    nn_dists  = np.full(N, np.inf)

    for i in range(N):
        lo        = max(0, i - window)
        hi        = min(N, i + window + 1)
        neighbors = np.concatenate([sorted_pts[lo:i], sorted_pts[i+1:hi]])
        if neighbors.size:
            nn_dists[i] = np.linalg.norm(neighbors - sorted_pts[i], axis=1).min()

    valid = nn_dists[np.isfinite(nn_dists) & (nn_dists > 1e-12)]
    if valid.size == 0:
        return 1.0

    d_min = valid.min()
    d_max = valid.max()

    if d_max < 1e-12:
        return 1.0

    return float(np.clip(d_min / d_max, 1e-3, 1.0))


# ------------------------------------------------------------------------------
# Argument parser (exposed so linspect wrapper can display subcommand help)
# ------------------------------------------------------------------------------

def get_parser():
    parser = argparse.ArgumentParser(
        prog='linspect bz',
        description='Plot the Brillouin-zone k-mesh from a LinReTraCe HDF5 file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
  linspect bz run.hdf5                       # kx-ky plane, auto scale
  linspect bz run.hdf5 --voronoi                  # with BZ boundary
  linspect bz run.hdf5 --axes x z --voronoi      # kx-kz plane with BZ boundary
  linspect bz run.hdf5 --scale 0.5          # manual marker scale
  linspect bz run.hdf5 --out bz.png
""")

    parser.add_argument('hdf5file',
                        help='Path to the LinReTraCe HDF5 input file.')
    parser.add_argument('--axes', nargs=2, default=['x', 'y'],
                        metavar=('AX1', 'AX2'),
                        choices=['x', 'y', 'z'],
                        help='Two reciprocal-space axes to plot: x, y, or z. Default: x y')
    parser.add_argument('--voronoi', action='store_true', default=False,
                        help='Draw the first Brillouin zone boundary as the Voronoi/Wigner-Seitz '
                             'cell of the reciprocal lattice. Off by default.')
    parser.add_argument('--out', default=None, metavar='FILE',
                        help='Save figure to FILE (e.g. bz.png). '
                             'If omitted, an interactive window is shown.')
    parser.add_argument('--cmap', default='viridis', metavar='NAME',
                        help='Matplotlib colormap for the out-of-plane colour coding. '
                             'Default: viridis')
    parser.add_argument('--scale', type=float, default=None, metavar='SCALE',
                        help='Marker size scaling factor in (0, 1]. '
                             '1 = standard size (120 pt^2), smaller values shrink symbols. '
                             'If omitted, the scale is estimated automatically from the '
                             'minimum k-point spacing (O(N log N), works for refined meshes).')
    return parser


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main(argv=None):
    parser = get_parser()
    args   = parser.parse_args(argv)

    if args.scale is not None and not (0 < args.scale <= 1.0):
        parser.error('--scale must be in the range (0, 1].')

    ax1_name, ax2_name = args.axes[0], args.axes[1]
    if ax1_name == ax2_name:
        parser.error('--axes must name two different directions.')

    ax1 = AXIS_INDEX[ax1_name]
    ax2 = AXIS_INDEX[ax2_name]
    ax3 = ({0, 1, 2} - {ax1, ax2}).pop()
    ax3_name = {v: k for k, v in AXIS_INDEX.items()}[ax3]

    # -- Read HDF5 -------------------------------------------------------------
    with h5py.File(args.hdf5file, 'r') as f:
        kvec      = f['/.unitcell/kvec'][:]    # (3,3) rows = b1, b2, b3
        kpts_frac = f['/.kmesh/points'][:]     # (N,3) fractional

    # -- Convert to Cartesian --------------------------------------------------
    kpts_cart  = frac_to_cart(kpts_frac, kvec)
    horiz      = kpts_cart[:, ax1]
    vert       = kpts_cart[:, ax2]
    color_vals = kpts_cart[:, ax3]

    # -- Reciprocal unit-cell parallelogram ------------------------------------
    b1_2d = np.array([kvec[0, ax1], kvec[0, ax2]])
    b2_2d = np.array([kvec[1, ax1], kvec[1, ax2]])
    para  = np.array([[0, 0], b1_2d, b1_2d + b2_2d, b2_2d, [0, 0]])

    # -- BZ boundary (optional) ------------------------------------------------
    if args.voronoi:
        bz_verts = bz_polygon_2d(kvec, ax1, ax2)

    # -- Marker size -----------------------------------------------------------
    s_max = 120   # pt^2 at scale = 1 on a uniform mesh

    if args.scale is not None:
        # Manual override
        scale = args.scale
        scale_source = f'manual (--scale {scale})'
    else:
        # Automatic: O(N log N) from sorted consecutive distances
        scale = auto_scale(kpts_cart)
        scale_source = f'auto ({scale:.3f})'

    marker_size = s_max * scale ** 2

    # -- Plot ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6, 6))

    # Reciprocal unit cell
    ax.plot(para[:, 0], para[:, 1],
            color='red', lw=1.85, ls='--', zorder=1, label='Recip. unit cell')

    # BZ polygon
    if args.voronoi:
        ax.fill(bz_verts[:-1, 0], bz_verts[:-1, 1],
                facecolor='#e8f0fb', edgecolor='#2a6abb', lw=2.0, zorder=2)
        ax.plot(bz_verts[:, 0], bz_verts[:, 1],
                color='#2a6abb', lw=2.0, zorder=3, label='1st Brillouin zone (Voronoi)')

    # k-points
    sc = ax.scatter(horiz, vert, c=color_vals, cmap=args.cmap,
                    s=marker_size, zorder=4, edgecolors='k', linewidths=0.6,
                    label=f'k-points (colour = k{ax3_name})')

    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(f'{AXIS_LABELS[ax3_name]}  (A\u207b\u00b9)', fontsize=11)

    # Reciprocal-vector arrows
    arrow_kw = dict(arrowstyle='->', color='black', lw=1.2, mutation_scale=12, zorder=5)
    for vec_idx, (dx, dy) in enumerate([b1_2d, b2_2d]):
        ax.annotate('', xy=(dx, dy), xytext=(0, 0), arrowprops=arrow_kw)
        ax.text(dx * 1.08, dy * 1.08,
                rf'$\mathbf{{b}}_{vec_idx + 1}$',
                ha='center', va='center', fontsize=11, color='#333333')

    # -- Axis limits -----------------------------------------------------------
    ax.set_aspect('equal')
    ax.autoscale(enable=False)

    all_x = list(horiz) + [b1_2d[0], b2_2d[0], 0.0]
    all_y = list(vert)  + [b1_2d[1], b2_2d[1], 0.0]
    if args.voronoi:
        all_x += list(bz_verts[:, 0])
        all_y += list(bz_verts[:, 1])

    xlo, xhi = min(all_x), max(all_x)
    ylo, yhi = min(all_y), max(all_y)
    xpad = 0.20 * (xhi - xlo)
    ypad = 0.20 * (yhi - ylo)
    ax.set_xlim(xlo - xpad, xhi + xpad)
    ax.set_ylim(ylo - ypad, yhi + ypad)

    # -- Legend placement ------------------------------------------------------
    # With --voronoi there are 3 legend entries and the BZ polygon fills much
    # of the axes area, making any in-axes corner placement prone to overlap.
    # The reliable solution is to place the legend outside the axes (below),
    # which also frees the full plot area for the data.
    # Without --voronoi (2 entries, no polygon) the standard in-axes placement
    # based on the b1+b2 resultant direction works fine.
    legend_kw: dict = dict(fontsize=9, framealpha=0.85)
    if args.voronoi:
        legend_kw.update(
            loc='upper center',
            bbox_to_anchor=(0.5, -0.12),
            bbox_transform=ax.transAxes,
            ncol=2,
            borderaxespad=0,
        )
    else:
        mean_vec = b1_2d + b2_2d
        angle = np.degrees(np.arctan2(mean_vec[1], mean_vec[0]))
        if -90 <= angle < 90:
            legend_loc = 'lower left' if angle >= 0 else 'upper left'
        else:
            legend_loc = 'upper right' if -180 <= angle < -90 else 'lower right'
        legend_kw['loc'] = legend_loc

    # -- Cosmetics -------------------------------------------------------------
    ax.axhline(0, color='#cccccc', lw=0.7, zorder=0)
    ax.axvline(0, color='#cccccc', lw=0.7, zorder=0)
    ax.set_xlabel(f'{AXIS_LABELS[ax1_name]}  (A\u207b\u00b9)', fontsize=13)
    ax.set_ylabel(f'{AXIS_LABELS[ax2_name]}  (A\u207b\u00b9)', fontsize=13)
    ax.set_title(
        f'Brillouin zone  \u2014  k{ax1_name}/k{ax2_name} plane\n'
        f'{len(kpts_frac)} k-points  |  marker scale: {scale_source}',
        fontsize=11)
    ax.legend(**legend_kw)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.tick_params(which='both', direction='in', top=True, right=True)

    plt.tight_layout()

    if args.out:
        fig.savefig(args.out, dpi=180, bbox_inches='tight')
        print(f'Figure saved to {args.out}')
    else:
        plt.show()


if __name__ == '__main__':
    main()
