#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
structure/checkfit.py : Visual postprocessing diagnostic for the
Sigma(0+)/Gamma/Z extraction performed by structure/w2dyn.py.

Plots the raw positive-Matsubara-frequency self-energy Sigma_LL'(i wm)
(real and imaginary parts) against the polynomial fit *actually used*
by w2dyn.extract_selfenergy -- via w2dyn.fit_diagnostics, which shares
its low-level fitting routine with the production code, so there is no
risk of the plot silently drifting from what gets written to the
scattering HDF5.

Designed as a postprocessing tool, run separately from the (potentially
expensive, cluster-batch) lwann --dmft extraction itself:
  * always saves figures to disk (PNG by default), so it works
    unattended/headless (no display required, e.g. over ssh or in a
    batch job);
  * --show additionally pops up an interactive window, if a display is
    available;
  * accepts MULTIPLE w2dynamics HDF5 files at once (e.g. a temperature
    scan) and overlays them in a single figure per orbital pair, which
    is exactly the comparison --fit-wmax is meant to make apples-to-
    apples (same energy window, different number of Matsubara points).
  * needs only the w2dynamics HDF5 file(s) -- no Wannier90 folder, no
    energy HDF5, no k-mesh -- so it can be run standalone, immediately
    after (or even during) a DMFT run, well before any lwann --dmft
    call.

Usage
-----
  python -m structure.checkfit FILE.hdf5 [FILE2.hdf5 ...] \\
      --elements diag \\
      --n-fit 20 --fit-degree 2 [--fit-wmax 0.3] \\
      --n-context 60 [--outdir ./fit_check] [--show]

  --elements accepts 'diag' (default, all N diagonal elements),
  'offdiag' (every unique off-diagonal pair, only meaningful for
  siw-full/offdiag runs), 'all' (diag + offdiag together), or one or
  more explicit 'L,Lp' pairs, e.g. --elements 0,0 0,1 1,1.
  For any off-diagonal pair (L != Lp), BOTH Sigma_LL'(i wm) and its
  partner Sigma_L'L(i wm) are fit independently and overlaid in the
  SAME figure (solid filled circles/line = (L,Lp), dashed hollow
  squares/line = (Lp,L)) -- their mutual agreement is the actual
  question of interest for an off-diagonal element (Sigma_LL' ~=
  Sigma_L'L is expected, NOT Hermitian conjugation, see the
  accompanying tex sec:offdiag), so this is more informative than
  plotting either direction alone.  Asking for just one of (L,Lp) or
  (Lp,L) (e.g. --elements 0,1) still produces this combined view.
  NOTE: --elements input is 0-indexed (plain Python/array convention,
  orbital 1 = '0'), but titles, axis labels, filenames, and console
  messages report 1-indexed orbitals (orbital 1 = "1"), the usual
  physics convention.  E.g. --elements 0,1 produces a file named
  sigma_fit_L1_Lp2.png.

Can also be imported and called directly, e.g. from lwann's
--check-fit shortcut:
  from structure.checkfit import make_plot
  make_plot(files=[...], elements_spec='diag', n_fit=10, fit_degree=2,
            fit_wmax=None, outdir='./fit_check', show=False)
"""

import argparse
import os
from typing import List, Optional, Sequence, Tuple, Union

import numpy as np

from structure import w2dyn


def _parse_elements(spec: Union[str, Sequence[str]], N: int) -> List[Tuple[int, int]]:
    """Turn 'diag' / 'offdiag' / 'all' / ['0,0', '0,1', ...] into a list of
    (L, Lp) tuples.  Input is 0-indexed (unchanged Python/array convention);
    only the *display* layer (titles, labels, filenames, console
    messages -- see make_plot below) reports 1-indexed orbitals.

    For 'offdiag' and 'all', only unique unordered off-diagonal pairs
    (L < Lp) are returned: make_plot always overlays the symmetric
    partner (Lp, L) automatically in the same figure whenever L != Lp
    (see the module docstring for why that overlay, not just the single
    element, is the diagnostically useful thing to look at for
    off-diagonal elements), so listing (L,Lp) and (Lp,L) as two separate
    requests here would just produce two duplicate figures.
    """
    if spec == 'diag' or spec is None:
        return [(L, L) for L in range(N)]
    if spec == 'offdiag':
        return [(L, Lp) for L in range(N) for Lp in range(L + 1, N)]
    if spec == 'all':
        return ([(L, L) for L in range(N)]
                + [(L, Lp) for L in range(N) for Lp in range(L + 1, N)])
    pairs = []
    for tok in spec:
        try:
            L_str, Lp_str = tok.split(',')
            pairs.append((int(L_str), int(Lp_str)))
        except ValueError:
            raise ValueError(
                "--elements entries must be 'diag', 'all', or 'L,Lp' "
                "(e.g. '0,1'); got {!r}".format(tok))
    return pairs


def make_plot(files: Sequence[str],
              elements_spec: Union[str, Sequence[str]] = 'diag',
              n_fit: int = 10,
              fit_degree: int = 2,
              fit_wmax: Optional[float] = None,
              spin_average: Optional[bool] = None,
              n_context: int = 60,
              y_padding: Optional[float] = 0.15,
              outdir: str = './fit_check',
              show: bool = False) -> List[str]:
    """
    Produce one figure per requested orbital pair, overlaying every input
    file (e.g. one figure per (L,Lp), one curve per temperature).

    Parameters mirror w2dyn.extract_selfenergy's (n_fit, fit_degree,
    fit_wmax, spin_average); see that function's docstring.  n_context
    sets how many low Matsubara points are shown beyond the fit window
    (in a lighter/hollow marker), for context.

    y_padding : float or None
        By default (0.15), each subplot's y-axis is clipped to the range
        spanned by the *raw data and the extracted intercept* shown,
        padded by this fraction of that range -- NOT by the full
        (possibly huge) excursion of the fitted polynomial curve outside
        its fit window, which is otherwise what matplotlib's automatic
        scaling would lock onto for high-degree fits and which makes it
        impossible to see how well the fit actually tracks the data near
        omega=0 (the only part that matters for Sigma0/Gamma/Z).  Use a
        larger value to zoom out a bit, or y_padding=None to disable
        clipping entirely and fall back to matplotlib's automatic range
        (i.e. the pre-existing behaviour, dominated by the fit curve).

    Returns the list of written figure paths.
    """
    import matplotlib
    if not show:
        matplotlib.use('Agg')   # headless-safe; must be set before pyplot import
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)

    datasets = [w2dyn.read_w2dyn(f) for f in files]
    N_list = [sum(d.orb_per_ineq) for d in datasets]
    if len(set(N_list)) != 1:
        raise ValueError(
            "Input files have different total orbital counts {}; "
            "they must be the same DMFT system at different "
            "temperatures.".format(dict(zip(files, N_list))))
    N = N_list[0]

    elements = _parse_elements(elements_spec, N)
    # normalise explicit user-typed pairs to canonical (L<=Lp) order and
    # drop duplicates, preserving first-seen order -- (1,0) and (0,1) are
    # the same request (both directions are always shown together below).
    elements = list(dict.fromkeys(
        (L, Lp) if L <= Lp else (Lp, L) for (L, Lp) in elements))

    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    written = []
    for (L, Lp) in elements:
        offdiag_pair = (L != Lp)
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
        any_plotted = False
        re_values, im_values = [], []   # raw data + star, for y-axis clipping

        for ifile, (d, f) in enumerate(zip(datasets, files)):
            if not d.offdiag and offdiag_pair:
                # off-diagonal element requested but this run only has
                # diagonal siw -- skip silently for this file, it simply
                # has no such data, rather than erroring out the whole
                # multi-file comparison.
                continue

            # One colour per FILE (temperature), explicit and fixed --
            # reused identically for every direction plotted below, so
            # colour always means "which file", never "which direction".
            color = color_cycle[ifile % len(color_cycle)]

            # Directions to overlay for this file: just (L,Lp) on the
            # diagonal; both (L,Lp) and (Lp,L) off-diagonal, since their
            # mutual agreement (Eq. offdiag_symmetry: Sigma_LL' ~= Sigma_L'L,
            # NOT Hermitian conjugation) is the actual question of interest
            # for an off-diagonal element -- a single direction alone
            # cannot show whether siw-full is converged/symmetric enough
            # for the symmetrisation step of extract_selfenergy to be
            # trustworthy.  Distinguished from each other by marker/
            # linestyle only (solid filled circle vs. dashed hollow square),
            # never by colour.
            directions = [(L, Lp)] if not offdiag_pair else [(L, Lp), (Lp, L)]

            for idir, (Ld, Lpd) in enumerate(directions):
                try:
                    info = w2dyn.fit_diagnostics(
                        d, Ld, Lpd, n_fit=n_fit, fit_degree=fit_degree,
                        fit_wmax=fit_wmax, spin_average=spin_average,
                        n_context=n_context)
                except ValueError as e:
                    print("  skipping {} for orbital pair ({},{}): {}".format(
                        f, Ld + 1, Lpd + 1, e))
                    continue
                any_plotted = True

                marker, mfc = ('o', color) if idir == 0 else ('s', 'none')
                ls = '-' if idir == 0 else '--'
                label = "{}  (T={:.1f} K), ({},{})".format(
                    os.path.basename(f), d.temperature, Ld + 1, Lpd + 1) \
                    if offdiag_pair else \
                    "{}  (T={:.1f} K)".format(os.path.basename(f), d.temperature)
                n_used = info['n_fit']

                for ax, comp, values in zip(axes, ('re', 'im'), (re_values, im_values)):
                    raw = info['siw'].real if comp == 're' else info['siw'].imag
                    fit = info['fit_re']   if comp == 're' else info['fit_im']
                    ax.plot(info['wm'][:n_used], raw[:n_used],
                            marker=marker, ls='None', ms=4,
                            mfc=mfc, mec=color, color=color,
                            label=label if comp == 're' else None)
                    values.extend(raw[:n_used])
                    # context points beyond the fit window: hollow, faded
                    if n_used < len(info['wm']):
                        ax.plot(info['wm'][n_used:], raw[n_used:],
                                marker=marker, ls='None', ms=4,
                                mfc='none', mec=color, alpha=0.4)
                        values.extend(raw[n_used:])
                    # fit curve, extrapolated down to wm=0 -- note this is
                    # NOT included in re_values/im_values: the y-range is
                    # clipped to the data (see y_padding above), not to
                    # whatever the polynomial does outside its fit window.
                    ax.plot(info['fit_wm'], fit, ls=ls, color=color, lw=1.5)
                    # cutoff used for the fit
                    ax.axvline(info['wm'][n_used - 1], color=color,
                               ls=':', lw=0.8, alpha=0.5)

                axes[0].plot(0.0,  info['Sigma0'], marker='*', ms=13,
                            color=color, markeredgecolor='k',
                            markeredgewidth=0.5, zorder=5)
                axes[1].plot(0.0, -info['Gamma'],  marker='*', ms=13,
                            color=color, markeredgecolor='k',
                            markeredgewidth=0.5, zorder=5)
                re_values.append(info['Sigma0'])
                im_values.append(-info['Gamma'])

        if not any_plotted:
            plt.close(fig)
            continue

        # Clip the y-axis to the raw-data (+ extracted intercept) range,
        # not to whatever the fitted polynomial does outside its fit
        # window -- see y_padding in the docstring above for why.
        if y_padding is not None:
            for ax, values in zip(axes, (re_values, im_values)):
                if not values:
                    continue
                vmin, vmax = min(values), max(values)
                span = vmax - vmin
                pad = span * y_padding if span > 0 else \
                    max(abs(vmin), abs(vmax), 1.0) * y_padding
                ax.set_ylim(vmin - pad, vmax + pad)

        # Display layer only, 1-based (physics convention): internal data
        # access above (fit_diagnostics(d, L, Lp, ...)) stays 0-indexed;
        # only what is shown/written here is shifted by +1.
        L1, Lp1 = L + 1, Lp + 1
        axes[0].set_xlabel(r'$\omega_m$ [eV]')
        axes[1].set_xlabel(r'$\omega_m$ [eV]')
        if offdiag_pair:
            axes[0].set_ylabel(r"Re $\Sigma_{LL'}(i\omega_m)$ [eV]")
            axes[1].set_ylabel(r"Im $\Sigma_{LL'}(i\omega_m)$ [eV]")
        else:
            axes[0].set_ylabel(r"Re $\Sigma_{%d%d}(i\omega_m)$ [eV]" % (L1, Lp1))
            axes[1].set_ylabel(r"Im $\Sigma_{%d%d}(i\omega_m)$ [eV]" % (L1, Lp1))
        axes[0].axhline(0, color='gray', lw=0.5)
        axes[1].axhline(0, color='gray', lw=0.5)
        axes[0].legend(fontsize=7.5, loc='best')
        if offdiag_pair:
            title = ("orbital pair ({},{})/({},{})  -- symmetry check "
                     "(solid filled=({},{}), dashed hollow=({},{}))   "
                     "fit_degree={}").format(
                L1, Lp1, Lp1, L1, L1, Lp1, Lp1, L1, fit_degree)
        else:
            title = 'orbital {}   fit_degree={}'.format(L1, fit_degree)
        title += '  fit_wmax={:.3f} eV'.format(fit_wmax) if fit_wmax is not None \
                 else '  n_fit={}'.format(n_fit)
        fig.suptitle(title, fontsize=10)
        fig.tight_layout(rect=[0, 0, 1, 0.94])

        outpath = os.path.join(
            outdir, 'sigma_fit_L{}_Lp{}.png'.format(L1, Lp1))
        fig.savefig(outpath, dpi=150)
        written.append(outpath)
        print("wrote", outpath)

        if show:
            plt.show()
        else:
            plt.close(fig)

    return written


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog='python -m structure.checkfit',
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('files', nargs='+', metavar='HDF5',
                    help='one or more w2dynamics HDF5 output files')
    p.add_argument('--elements', nargs='+', default=['diag'],
                    help="'diag' (default), 'offdiag' (every unique "
                         "off-diagonal pair), 'all' (both), or explicit "
                         "'L,Lp' pairs, e.g. --elements 0,0 0,1 1,1. "
                         "Off-diagonal pairs are always shown together "
                         "with their (Lp,L) partner in one figure (see "
                         "module docstring). 0-indexed input (orbital 1 "
                         "= '0'); output (titles/filenames/console) is "
                         "1-indexed.")
    p.add_argument('--n-fit', type=int, default=10,
                    help='number of lowest Matsubara points for the fit '
                         '(default: 10); ignored if --fit-wmax is given')
    p.add_argument('--fit-degree', type=int, default=2,
                    help='polynomial degree, >= 1 (default: 2, quadratic). '
                         'Higher orders can help for runs with very good '
                         'statistics, but need a correspondingly larger '
                         '--n-fit/--fit-wmax to stay well-conditioned -- '
                         'always check the figure (this tool) before '
                         'trusting a higher-order extrapolation.')
    p.add_argument('--fit-wmax', type=float, default=None, metavar='ENERGY',
                    help='[eV] use all Matsubara points up to this energy '
                         'cutoff instead of a fixed --n-fit count; '
                         'recommended when comparing several temperatures '
                         '(same as the lwann --dmft flag of the same name)')
    p.add_argument('--no-spin-average', action='store_true', default=False,
                    help='use the spin-0 channel only instead of '
                         'auto-averaging for paramagnetic runs')
    p.add_argument('--n-context', type=int, default=60,
                    help='number of low Matsubara points to show beyond '
                         'the fit window, for context (default: 60)')
    p.add_argument('--y-padding', type=float, default=0.15,
                    help='y-axis is clipped to the range of the raw data '
                         'shown (+ the extracted Sigma0/Gamma intercept), '
                         'padded by this fraction of that range (default: '
                         '0.15) -- NOT to the fitted curve, which can '
                         'range far more widely outside its fit window, '
                         'especially for higher fit_degree; use '
                         '--full-yrange to see that instead.')
    p.add_argument('--full-yrange', action='store_true', default=False,
                    help='disable y-axis clipping; let matplotlib '
                         'auto-scale to the full fitted curve as well '
                         '(the pre-existing behaviour). Useful to '
                         'deliberately inspect how badly a high-degree '
                         'fit misbehaves outside its window.')
    p.add_argument('--outdir', default='./fit_check',
                    help='output directory for the figures '
                         '(default: ./fit_check)')
    p.add_argument('--show', action='store_true', default=False,
                    help='also open an interactive window '
                         '(requires a display)')
    return p


def main(argv=None):
    args = _build_parser().parse_args(argv)
    elements_spec = (args.elements[0]
                      if len(args.elements) == 1 and args.elements[0] in ('diag', 'offdiag', 'all')
                      else args.elements)
    spin_average = False if args.no_spin_average else None
    y_padding = None if args.full_yrange else args.y_padding
    make_plot(
        files        = args.files,
        elements_spec= elements_spec,
        n_fit        = args.n_fit,
        fit_degree   = args.fit_degree,
        fit_wmax     = args.fit_wmax,
        spin_average = spin_average,
        n_context    = args.n_context,
        y_padding    = y_padding,
        outdir       = args.outdir,
        show         = args.show,
    )


if __name__ == '__main__':
    main()
