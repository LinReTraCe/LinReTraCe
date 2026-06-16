#! /usr/bin/env python3
"""
scripts/ltb/run.py
==================
``ltb run`` subcommand: build a LinReTraCe HDF5 file from a tight-binding
model on a regular k-mesh.

New in this version
-------------------
``--write-hr STEM`` / ``--write_hr STEM``
    After the HDF5 file has been written, export the TightBinding object as a
    Wannier90-compatible file set (``<STEM>_hr.dat``, ``<STEM>.nnkp``,
    ``<STEM>.wout``) so that the ``lwann`` input generator can read the same
    model without re-running ``ltb``.
"""
from __future__ import annotations

import sys
from pathlib import Path

# Ensure the package root is on sys.path when invoked directly by path.
# scripts/ltb/ sits two levels below the package root.
_root = Path(__file__).resolve().parents[2]
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

# ── OpenMP / BLAS thread count ────────────────────────────────────────────────
# Must happen before numpy (or any BLAS-linked library) is imported.
from scripts.ltb.openmp_utils import (
    preparse_openmp, add_openmp_argument, openmp_info_line,
)
_openmp_ncores = preparse_openmp()
# ─────────────────────────────────────────────────────────────────────────────

import os
import logging
import argparse
from argparse import RawTextHelpFormatter

import numpy as np

from structure.inout          import h5output
from structure.tb             import TightBinding
from structure.auxiliary      import LogFormatter
from structure.wannier_export import write_wannier90


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description='Build a LinReTraCe HDF5 file from a tight-binding model'
                    ' on a regular k-mesh.',
        formatter_class=RawTextHelpFormatter,
        prog='ltb run',
        epilog="""\
tb_file:
  begin hopping:      hopping parameter [eV]
                      sign convention e(k) ~ -t * (1-2.delta_(R,0).delta(l,l')) e^{ikR}
  begin atoms:        fractional atomic positions inside unit cell
                      # used to determine the unit cell symmetries only
  begin real_lattice: lattice vectors [Angstrom]
                      each row = 1 vector, each column = Cartesian x/y/z component
  optional:
  begin orbitals:     fractional atomic positions of the orbitals
                      referenced in 'begin hopping'.
                      # used for the multi-cell peierls correction

begin hopping
# x  y  z  orbital1 orbital2 hopping.real [hopping.imag]
  0 +1  0  1        1        0.25
  0 -1  0  1        1        0.25
 +1 +1  0  1        2        0.05 0.005
...
end hopping

begin atoms
# ineq (starting at 1) x y z
  1                    0 0 0
...
end atoms

begin real_lattice
# x y z
  5 0 0
  0 5 0
  0 0 1
end real_lattice

OPTIONAL:
begin orbitals
# orbital (from hopping) x y z
  1                      0 0 0
...
end orbitals
""")

    # ── mandatory ─────────────────────────────────────────────────────────────
    parser.add_argument('tb_file',  type=str,   help='tight binding file')
    parser.add_argument('nkx',      type=int,   help='number of k-points in x-direction')
    parser.add_argument('nky',      type=int,   help='number of k-points in y-direction')
    parser.add_argument('nkz',      type=int,   help='number of k-points in z-direction')
    parser.add_argument('filling',  type=float, help='number of electrons in the system'
                                                     ' (1 filled band = 2 electrons)')

    # ── optional ──────────────────────────────────────────────────────────────
    parser.add_argument('-o', '--output',   default=None,
                        help='Output file name')
    parser.add_argument('-p', '--plot',     default=False, action='store_true',
                        help='Plot bands / DOS / NOS')
    parser.add_argument('--kshift',         default=False, action='store_true',
                        help='shift k-grid by half a k-point')
    parser.add_argument('--mushift',        default=False, action='store_true',
                        help='shift energies such that the chemical potential is at mu = 0')
    parser.add_argument('--intraonly',      default=False, action='store_true',
                        help='only save intra-band elements')
    parser.add_argument('--intra',          type=float,
                        help='set all intra-band elements to given value')
    parser.add_argument('--inter',          type=float,
                        help='set all inter-band elements to given value')
    parser.add_argument('--red',            default=False, action='store_true',
                        help='make k-grid reducible')
    parser.add_argument('--mu',             type=float,
                        help='use provided chemical potential instead of provided filling'
                             ' (debugging purposes)')
    parser.add_argument('--corronly',       default=False, action='store_true',
                        help='only use the multi-orbital Peierls correction term'
                             ' (derivative term is set to 0)')
    parser.add_argument('--sparse_rotation', default=False, action='store_true',
                        help='Use sparse matrix multiplication when rotating velocities\n'
                             'and curvatures into the Kohn-Sham basis.\n'
                             'Speedup grows as ~N_orb/z for large supercells (z = coordination).\n'
                             'No effect on results. Recommended for N_orb > ~100.')
    parser.add_argument('--vector',         default=False, action='store_true',
                        help='save also the full hamiltonian and transformations to HDF5')

    # ── Wannier90 export ──────────────────────────────────────────────────────
    parser.add_argument(
        '--write-hr', '--write_hr',
        metavar='STEM',
        default=None,
        dest='write_hr',
        help=(
            'After writing the HDF5 file, export the TB model in\n'
            'Wannier90 format so that lwann can read it directly.\n'
            'Three files are created:\n'
            '  <STEM>_hr.dat  real-space Hamiltonian H(R)\n'
            '  <STEM>.nnkp    lattice vectors, k-mesh, projections\n'
            '  <STEM>.wout    length-unit stub\n'
            'The directory containing these files can then be passed\n'
            'to lwann as its system argument.  If STEM contains a\n'
            'directory component it will be created if necessary.\n'
            'Example:\n'
            '  ltb run model.tb 10 10 10 2.0 --red --write-hr wann\n'
            '  lwann . --charge 2.0\n'
            'Or with a subdirectory:\n'
            '  ltb run model.tb 10 10 10 2.0 --red --write-hr out/wann\n'
            '  lwann out/ --charge 2.0\n'
            'A round-trip test (dispersion + file format) is available at\n'
            '  tests/test_wannier_export.py'
        ),
    )

    parser.add_argument('--debug', help=argparse.SUPPRESS,
                        default=False, action='store_true')
    add_openmp_argument(parser)
    return parser.parse_args(args)


def main(argv=None):
    error = lambda msg: sys.exit('\nltb: {}'.format(msg))
    args  = parse_args(argv)
    debug = args.debug

    # ── logging ───────────────────────────────────────────────────────────────
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG if debug else logging.INFO)
    console = logging.StreamHandler()
    console.setFormatter(LogFormatter())
    console.setLevel(logging.DEBUG if debug else logging.INFO)
    logger.addHandler(console)

    logger.info(openmp_info_line(_openmp_ncores))

    # ── build TightBinding object ─────────────────────────────────────────────
    irr = not args.red
    try:
        tb = TightBinding(
            nkx=args.nkx,
            nky=args.nky,
            nkz=args.nkz,
            irreducible=irr,
            kshift=args.kshift,
        )
    except Exception as e:
        error(str(e) + '\nExit.')

    # ── compute dispersion, velocities, etc. ──────────────────────────────────
    try:
        tb.computeData(
            tbfile=args.tb_file,
            charge=args.filling,
            mu=args.mu,
            mushift=args.mushift,
            corronly=args.corronly,
            vector=args.vector,
            sparse_rotation=args.sparse_rotation,
        )
    except Exception as e:
        error(str(e) + '\nExit.')

    if args.intra is not None:
        tb.setDiagonal(args.intra)
    if args.inter is not None:
        tb.setOffDiagonal(args.inter)

    if args.intraonly:
        tb.bopticfull = False
        tb.opticfull  = False

    # ── HDF5 output ───────────────────────────────────────────────────────────
    fname = (args.output if args.output is not None
             else 'lrtc-tb-{}-{}-{}-{}.hdf5'.format(
                 tb.nkx, tb.nky, tb.nkz, 'irr' if irr else 'red'))
    h5output(fname, tb, tb)

    logger.info("Output file {!r} successfully created.".format(fname))
    logger.info("File size: {} MB.".format(os.stat(fname).st_size / 1_000_000))

    # ── Wannier90 export (optional) ───────────────────────────────────────────
    if args.write_hr is not None:
        stem = args.write_hr
        # Create parent directory if the stem contains a path component.
        parent = Path(stem).parent
        if str(parent) not in ('', '.'):
            parent.mkdir(parents=True, exist_ok=True)
        try:
            write_wannier90(tb, stem)
        except Exception as e:
            error('Wannier90 export failed: ' + str(e) + '\nExit.')

    # ── optional plot ─────────────────────────────────────────────────────────
    if args.plot:
        try:
            logging.getLogger('matplotlib').setLevel(logging.WARNING)
            import matplotlib.pyplot as plt
        except ImportError:
            error('--plot requires the matplotlib library')

        for iband in range(tb.energyBandMax):
            plt.plot(tb.energies[0][:, iband],
                     label='band {}'.format(iband + 1), lw=2)
        plt.axhline(tb.mu,
                    label='mu_tb = {:.3f}'.format(tb.mu),
                    color='black', lw=1, ls='--')
        plt.xlabel(r'$k_i$')
        plt.ylabel(r'$\varepsilon(k_i)$')
        plt.legend(loc='best')
        plt.show()

        tb.calcDOS(gamma=0.03, npoints=10_000, windowsize=1.5)

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twinx()
        for ispin in range(tb.spins):
            ax1.plot(tb.dosaxis, tb.dos[ispin],
                     label='dos', color='blue', lw=2)
            ax2.plot(tb.dosaxis, tb.nos[ispin],
                     label='nos', color='red',  lw=2)

        ax1.axvline(x=tb.mu,      color='black', lw=1, ls='-')
        ax1.set_ylim(ymin=0)
        ax1.set_ylabel(r'$\mathrm{dos}$')
        ax1.set_xlabel(r'$\mu$ [eV]')
        ax1.legend(loc='center left')

        ax2.axhline(y=tb.charge,  color='black', lw=1, ls='-')
        ax2.set_ylim(ymin=0)
        ax2.set_ylabel(r'$\mathrm{nos}$')
        ax2.legend(loc='center right')

        plt.show()


if __name__ == '__main__':
    main()
