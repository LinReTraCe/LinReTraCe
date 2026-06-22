#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
structure/w2dyn.py
==================
Reader and self-energy fitting module for w2dynamics DMFT output files.

This module is the only part of the lwann/LinReTraCe interface that knows
about the w2dynamics HDF5 format.  Everything it exports is expressed in
plain NumPy arrays so that the calling code in wannier_dmft.py remains
solver-agnostic.

Public API
----------
  read_w2dyn(path)              -> W2dynData
  extract_selfenergy(data, ...) -> SelfEnergyParams
  read_and_sort_temperatures([paths], ...) -> [SelfEnergyParams]

Self-energy decomposition
--------------------------
The DMFT self-energy in the orbital basis is linearised as

    Sigma_LL'(iw_m) ~ Sigma0_LL'  -  i*Gamma_LL'  +  i*w_m*(delta_LL' - Zinv_LL')

where Sigma0_LL', Gamma_LL', and Zinv_LL' are all frequency-independent
coefficients extracted by polynomial fit to the lowest positive Matsubara
frequencies.

  Sigma0_orb : Sigma0_LL' = Re[Sigma_LL'(0+)]   -- REAL for every (L,L'),
               on AND off the diagonal.  (Stored as complex128 for
               interface uniformity with H(k); the imaginary part is
               always exactly 0 by construction, never populated.)
  Gamma_orb  : Gamma_LL' = -Im[Sigma_LL'(0+)]   -- REAL for every (L,L'),
               on AND off the diagonal.  Gamma_LL >= 0 is the physical
               causality requirement; the full matrix is expected to be
               (to the same Fermi-liquid/causality accuracy used
               throughout this module) SYMMETRIC and positive
               semi-definite, NOT antisymmetric.
  Zinv_orb   : DIAGONAL (from diagonal Fermi-liquid slope fit only)

It is *not* correct to keep the full complex Sigma_LL'(0+) (Re+iIm) in
Sigma0_orb for L!=L' and restrict Gamma_orb to the diagonal.  The naive
argument for that ("Sigma_LL'(iwn)^* = Sigma_L'L(-iwn) implies Sigma0_orb
is Hermitian") is wrong: it silently extends the POSITIVE-frequency fit
of element (L',L) to NEGATIVE frequencies, which is invalid whenever
that element has a finite Gamma_L'L (the one-sided iwn->0+ and iwn->0-
limits of Im[Sigma_L'L] generically differ by 2*Gamma_L'L; that is
exactly the standard Im[Sigma]->-Gamma*sgn(wn)-type jump at wn=0).  The
two independently-fitted intercepts therefore satisfy, to the accuracy
of the present Fermi-liquid extrapolation,
    Sigma_LL'(0+) ~= Sigma_L'L(0+)        (plain SYMMETRY, no conjugation)
not the Hermitian relation Sigma_LL'(0+) ~= conj(Sigma_L'L(0+)).  Since
Gamma_LL' is, by its very definition (-Im of a single complex number),
already a REAL matrix, "Hermitian" and "symmetric" coincide for it, and
the standard causality/positive-semi-definite broadening-matrix
argument that already justifies Gamma_LL>=0 on the diagonal generalises
directly to: the full matrix Gamma_orb should be symmetric, real,
PSD -- not antisymmetric (which is what would be required for i*Im(.)
restricted off-diagonal to itself be Hermitian).  Putting the
off-diagonal Im[Sigma(0+)] into H_tilde (with a "+i" sign, as if it
were Hermitian) therefore generically BREAKS the Hermiticity of
H_tilde.  The fix used here keeps the split clean and automatically
Hermitian-safe: Re[Sigma(0+)] (real symmetric, for every L,L') goes
into Sigma0_orb / H_tilde; -Im[Sigma(0+)] (real, for every L,L') goes
into Gamma_orb, which is then projected and truncated to its diagonal
in the QP basis exactly as already described for Approximation 3 (see
the accompanying tex, sec:QPdiag) -- no special-casing of the diagonal
vs. off-diagonal orbital case is needed any more.

Sign convention for double-counting (PROVISIONAL):
  G^{-1}(k,iw) = (iw+mu)*I - H(k) - Sigma(iw) - dc
  so dc is ADDED to Sigma0 when forming H_tilde.

w2dynamics storage conventions
--------------------------------
  - siw shape: (norb, nspins, 2*Niw): first Niw = negative freqs (descending),
    last Niw = positive freqs (ascending).  Only the positive half is used.
  - siw-full shape: (norb1, spin1, norb2, spin2, 2*Niw) for qmc.offdiag != 0.
  - Im Sigma(iw_n>0) < 0 in the standard causal convention.
  - Gamma_LL = -Im Sigma_LL(iw_n->0+) > 0.
  - Config metadata: .config.attrs stores general.magnetism, general.beta, etc.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import h5py

logger = logging.getLogger(__name__)

_KB_EV = 8.61733034e-5   # Boltzmann constant [eV / K]


# ─────────────────────────────────────────────────────────────────────────────
# Data containers
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class W2dynData:
    """
    Raw data extracted from one w2dynamics HDF5 output file.

    Attributes
    ----------
    path : str
    beta : float             Inverse temperature [eV^-1].
    temperature : float      Temperature [K].
    mu : float               DMFT chemical potential [eV].
    magnetism : str          'para', 'ferro', etc.
    offdiag : bool           True when qmc.offdiag != 0.
    n_ineq : int             Number of inequivalent correlated sites.
    orb_per_ineq : list[int] Number of active orbitals per site.
    pos_iw : ndarray (Niw,)  Positive Matsubara frequencies [eV].
    siw_ineq : list[ndarray]
        Per site.  offdiag=False: (norb, nspins, Niw) complex.
                   offdiag=True : (norb, nspins, norb, nspins, Niw) complex.
        Positive-frequency half already extracted.
    dc_ineq : list[ndarray]  Per site, shape (norb, nspins) [eV].
    smom_ineq : list[ndarray] Per site, shape (norb, nspins, 2).
    """
    path: str
    beta: float
    temperature: float
    mu: float
    magnetism: str
    offdiag: bool
    n_ineq: int
    orb_per_ineq: list
    pos_iw: np.ndarray
    siw_ineq: list
    dc_ineq: list
    smom_ineq: list


@dataclass
class SelfEnergyParams:
    """
    Linearised QP parameters extracted from one DMFT run.

    Attributes
    ----------
    temperature : float      [K]
    beta : float             [eV^-1]
    mu : float               DMFT chemical potential [eV]
    n_orb : int              Total Wannier orbital count.
    spin_average : bool
    Zinv_orb : ndarray (n_orb, n_orb) float
        Diagonal inverse QP weight.  Off-diagonal entries are 0.
    Sigma0_orb : ndarray (n_orb, n_orb) complex
        Real-symmetric matrix Re[Sigma_LL'(0+)], for every (L,L') pair.
        Stored as complex128 for interface uniformity with H(k); the
        imaginary part is always exactly 0 by construction.
    Gamma_orb : ndarray (n_orb, n_orb) float
        Scattering-rate matrix, Gamma_LL' = -Im[Sigma_LL'(0+)], for
        every (L,L') pair (on AND off the diagonal for offdiag=True
        runs; off-diagonal entries are 0 for offdiag=False runs, where
        siw alone carries no off-diagonal information).  Expected to be
        symmetric and positive semi-definite; truncation to its
        diagonal in the QP basis ("Approximation 3") happens later,
        downstream in wannier_dmft.py, not here.
    dc_orb : ndarray (n_orb,) float
        Spin-averaged double-counting per orbital [eV].
    n_fit : int
    fit_degree : int
    fit_wmax : float or None  If not None, energy cutoff used for n_fit.
    """
    temperature: float
    beta: float
    mu: float
    n_orb: int
    spin_average: bool
    Zinv_orb: np.ndarray
    Sigma0_orb: np.ndarray
    Gamma_orb: np.ndarray
    dc_orb: np.ndarray
    n_fit: int
    fit_degree: int
    fit_wmax: object  # float or None


# ─────────────────────────────────────────────────────────────────────────────
# Reader
# ─────────────────────────────────────────────────────────────────────────────

def read_w2dyn(path: str) -> W2dynData:
    """Read a w2dynamics HDF5 and return a W2dynData container."""
    path = str(path)
    logger.info("Reading w2dynamics HDF5: {}".format(path))

    with h5py.File(path, 'r') as h5:
        cfg = dict(h5['.config'].attrs)

        beta        = float(cfg['general.beta'])
        temperature = 1.0 / (beta * _KB_EV)

        # mu: read the actual converged chemical potential of the LAST
        # iteration, not the static config attribute.  general.mu (in
        # .config.attrs) is merely the value Parameters.in was started
        # with; whenever general.mu_search is active (dichotomy / nloc
        # density search) the converged mu drifts away from it over the
        # DMFT loop and only dmft-last/mu/value reflects what was
        # actually used to obtain the siw/dc/smom data we are about to
        # read.  (They happen to coincide for single-shot,
        # DMFTsteps=1 test runs, which can mask this bug.)
        mu_cfg = float(cfg['general.mu'])
        if 'dmft-last/mu/value' in h5:
            mu = float(h5['dmft-last/mu/value'][()])
            if abs(mu - mu_cfg) > 1e-6:
                logger.warning(
                    "general.mu (config) = {:.6f} eV differs from "
                    "dmft-last/mu/value = {:.6f} eV; using the converged "
                    "dmft-last value.".format(mu_cfg, mu))
        else:
            logger.warning(
                "dmft-last/mu/value not found; falling back to "
                "general.mu (config) = {:.6f} eV. This may be stale if "
                "a chemical-potential search was used.".format(mu_cfg))
            mu = mu_cfg

        magnetism   = str(cfg['general.magnetism'])
        offdiag     = (int(cfg['qmc.offdiag']) != 0)
        n_ineq      = int(cfg['general.nat'])

        orb_per_ineq = []
        for iat in range(1, n_ineq + 1):
            nd  = int(cfg.get('atoms.{}.nd'.format(iat), 0))
            np_ = int(cfg.get('atoms.{}.np'.format(iat), 0))
            orb_per_ineq.append(nd + np_)

        logger.info("  beta={:.4f} eV^-1  T={:.2f} K  mu={:.5f} eV  "
                    "magnetism={}  offdiag={}  n_ineq={}".format(
                        beta, temperature, mu, magnetism, offdiag, n_ineq))
        logger.info("  orbitals per ineq: {}".format(orb_per_ineq))

        # Positive Matsubara frequencies.  .axes/pos-iw (a precomputed
        # positive-only half) is written by some w2dynamics runs but not
        # all -- it is absent, e.g., from runs that also compute two-particle
        # quantities (where the .axes group instead holds iwf-p3, iwb-p2,
        # iwb-g4, etc., the vertex/susceptibility frequency axes, and no
        # pos-iw at all).  .axes/iw -- the full (negative+positive)
        # one-particle Matsubara axis -- is always present and unambiguous:
        # it is sorted ascending, symmetric about 0, with exactly as many
        # negative as positive entries, so its upper half is exactly the
        # positive-frequency axis that siw/siw-full are already sliced
        # against below.  Derive pos_iw from it directly; cross-check
        # against .axes/pos-iw when that also happens to be present.
        iw_full   = h5['.axes/iw'][()].astype(np.float64)
        Niw       = len(iw_full) // 2
        pos_iw    = iw_full[Niw:]
        if len(iw_full) % 2 != 0:
            raise ValueError(
                ".axes/iw has odd length {}; expected an even number of "
                "Matsubara points (symmetric about 0).".format(len(iw_full)))
        if '.axes/pos-iw' in h5:
            pos_iw_ref = h5['.axes/pos-iw'][()].astype(np.float64)
            if not np.allclose(pos_iw, pos_iw_ref, rtol=1e-10, atol=1e-12):
                logger.warning(
                    ".axes/pos-iw present but does not match the upper "
                    "half of .axes/iw; using .axes/iw (the one always "
                    "guaranteed to match siw's frequency slicing below).")

        siw_ineq, dc_ineq, smom_ineq = [], [], []

        for iat in range(1, n_ineq + 1):
            grp = h5['dmft-last/ineq-{:03d}'.format(iat)]

            dc   = grp['dc/value'][()].astype(np.float64)
            smom = grp['smom/value'][()].astype(np.float64)
            dc_ineq.append(dc)
            smom_ineq.append(smom)

            if offdiag:
                raw = grp['siw-full/value'][()].astype(np.complex128)
                siw_pos = raw[..., Niw:]          # positive half: last Niw entries
            else:
                raw = grp['siw/value'][()].astype(np.complex128)
                siw_pos = raw[:, :, Niw:]

            siw_ineq.append(siw_pos)
            logger.info("  ineq-{:03d}: dc shape {}  siw+ shape {}".format(
                iat, dc.shape, siw_pos.shape))

    return W2dynData(
        path=path, beta=beta, temperature=temperature, mu=mu,
        magnetism=magnetism, offdiag=offdiag, n_ineq=n_ineq,
        orb_per_ineq=orb_per_ineq, pos_iw=pos_iw,
        siw_ineq=siw_ineq, dc_ineq=dc_ineq, smom_ineq=smom_ineq,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Merge inequivalent sites → full Wannier basis
# ─────────────────────────────────────────────────────────────────────────────

def _merge_ineq(data: W2dynData) -> tuple:
    """Block-diagonal merge of per-site siw and dc arrays."""
    N_total = sum(data.orb_per_ineq)
    nspins  = data.siw_ineq[0].shape[1]
    Niw     = len(data.pos_iw)

    if data.offdiag:
        siw_full = np.zeros((N_total, nspins, N_total, nspins, Niw), dtype=np.complex128)
        dc_full  = np.zeros((N_total, nspins), dtype=np.float64)
        i0 = 0
        for iat, norb in enumerate(data.orb_per_ineq):
            i1 = i0 + norb
            siw_full[i0:i1, :, i0:i1, :, :] = data.siw_ineq[iat]
            dc_full[i0:i1, :]                = data.dc_ineq[iat]
            i0 = i1
    else:
        siw_full = np.zeros((N_total, nspins, Niw), dtype=np.complex128)
        dc_full  = np.zeros((N_total, nspins), dtype=np.float64)
        i0 = 0
        for iat, norb in enumerate(data.orb_per_ineq):
            i1 = i0 + norb
            siw_full[i0:i1, :, :] = data.siw_ineq[iat]
            dc_full[i0:i1, :]     = data.dc_ineq[iat]
            i0 = i1

    return siw_full, dc_full


# ─────────────────────────────────────────────────────────────────────────────
# Fitting utilities
# ─────────────────────────────────────────────────────────────────────────────

def _n_fit_from_wmax(wm: np.ndarray, fit_wmax: float, beta: Optional[float] = None) -> int:
    """Return the number of Matsubara points with wm <= fit_wmax."""
    n = int(np.searchsorted(wm, fit_wmax, side='right'))
    if n < 2:
        beta_str = "{:.1f}".format(beta) if beta is not None else "unknown"
        raise ValueError(
            "fit_wmax={} eV captures fewer than 2 Matsubara points "
            "(beta={}, wm[0]={:.4f}).  Increase fit_wmax or reduce beta."
            .format(fit_wmax, beta_str, wm[0]))
    return n


def _polyfit_element(siw_ll: np.ndarray, wm: np.ndarray,
                      n_fit: int, fit_degree: int) -> tuple:
    """
    Low-level least-squares polynomial fit of one Matsubara element
    Sigma_LL'(iw_n), shared by _fit_element (production extraction) and
    fit_diagnostics (visual inspection), so the diagnostic plot is
    guaranteed to show exactly the curve that gets used, never a
    re-implementation that could silently drift from it.

    Returns
    -------
    re_c, im_c : ndarray, ndarray
        numpy.polyfit coefficient arrays (highest degree first) for the
        real and imaginary parts respectively, fit to the lowest n_fit
        Matsubara points.
    """
    wm_lo = wm[:n_fit]
    im_c  = np.polyfit(wm_lo, siw_ll[:n_fit].imag, fit_degree)
    re_c  = np.polyfit(wm_lo, siw_ll[:n_fit].real, fit_degree)
    return re_c, im_c


def _fit_element(siw_ll: np.ndarray, wm: np.ndarray,
                 n_fit: int, fit_degree: int,
                 diagonal: bool) -> tuple:
    """
    Fit one Matsubara element Sigma_LL'(iw_n) at n_fit lowest frequencies.

    Returns
    -------
    Zinv_ll  : float   1 - d(Im Sigma)/d(wm) |_{wm->0}.
               Meaningful only for diagonal elements (L==L'); callers
               must ignore/overwrite this for L!=L' (Zinv_orb stays
               diagonal -- off-diagonal Z^-1 is not modelled here).
    Sigma0_ll: complex Re[Sigma_LL'(0+)] ONLY, imaginary part forced to
               exactly 0.  This holds for diagonal AND off-diagonal
               (L,L') pairs alike: see the module docstring for why
               keeping the off-diagonal Im[Sigma(0+)] in Sigma0
               instead would be incorrect -- it would generically
               break the Hermiticity of H_tilde.
    Gamma_ll : float   -Im[Sigma_LL'(0+)].  >= 0 expected on the
               diagonal (causality); off-diagonal entries are real
               numbers of either sign in principle, but the full
               Gamma_orb matrix assembled from these is expected to be
               symmetric and positive semi-definite (see module
               docstring) -- this is checked, not assumed, by the
               caller.
    """
    re_c, im_c = _polyfit_element(siw_ll, wm, n_fit, fit_degree)

    # np.polyfit returns coefficients highest-power-first: for degree p,
    # c[-1] is the constant (intercept) term and c[-2] is the linear
    # (slope) term, regardless of p -- this generalises cleanly to any
    # fit_degree >= 1, not just 1 or 2.
    if fit_degree < 1:
        raise ValueError("fit_degree must be >= 1, got {}".format(fit_degree))
    slope_im, intercept_im = im_c[-2], im_c[-1]
    intercept_re           = re_c[-1]

    Zinv_ll = 1.0 - slope_im if diagonal else 1.0

    # Sigma0_ll carries ONLY the real part, for every (L,L') pair,
    # diagonal or off-diagonal.  Im[Sigma_LL'(0+)] is returned
    # separately as Gamma_ll for every (L,L') pair and must never be
    # folded into Sigma0_ll: H_tilde = H(k) + Sigma0_orb + ... is only
    # guaranteed Hermitian if Sigma0_orb is real-symmetric.
    Sigma0_ll = complex(intercept_re, 0.0)
    Gamma_ll  = -intercept_im   # >=0 expected on the diagonal (causality)

    return float(Zinv_ll), Sigma0_ll, float(Gamma_ll)


# ─────────────────────────────────────────────────────────────────────────────
# Shared spin-averaging helper (used by extract_selfenergy and fit_diagnostics)
# ─────────────────────────────────────────────────────────────────────────────

def _get_siw_work(data: 'W2dynData', spin_average: Optional[bool]) -> tuple:
    """
    Resolve spin-averaging exactly as extract_selfenergy does, and return
    (siw_work, dc_work, N_total, spin_average_used).  Factored out so
    fit_diagnostics reuses precisely the same data that gets fitted in
    production, rather than re-deriving it.
    """
    if spin_average is None:
        spin_average = (data.magnetism.lower() == 'para')

    siw_full, dc_full = _merge_ineq(data)
    N_total = sum(data.orb_per_ineq)
    nspins  = siw_full.shape[1]

    if spin_average:
        if data.offdiag:
            siw_work = np.zeros_like(siw_full[:, 0, :, 0, :])
            for s in range(nspins):
                siw_work += siw_full[:, s, :, s, :]
            siw_work /= nspins
        else:
            siw_work = np.mean(siw_full, axis=1)
        dc_work = np.mean(dc_full, axis=1)
    else:
        if data.offdiag:
            siw_work = siw_full[:, 0, :, 0, :]
        else:
            siw_work = siw_full[:, 0, :]
        dc_work = dc_full[:, 0]

    return siw_work, dc_work, N_total, spin_average


def fit_diagnostics(data: 'W2dynData', L: int, Lp: int,
                     n_fit: int = 10, fit_degree: int = 2,
                     fit_wmax: Optional[float] = None,
                     spin_average: Optional[bool] = None,
                     n_context: int = 60, n_curve: int = 200) -> dict:
    """
    Return raw Matsubara data and the corresponding fit curve for one
    orbital pair (L, Lp), for visual inspection (used by
    structure/checkfit.py).  Calls the exact same low-level fitting
    routine (_polyfit_element) as extract_selfenergy, so what is plotted
    is guaranteed to be what gets used in the real DMFT extraction --
    there is no separate/duplicated fitting logic here.

    Parameters mirror extract_selfenergy's (n_fit, fit_degree, fit_wmax,
    spin_average); see that function's docstring.

    Returns a dict with keys:
      wm        : (n_context,) Matsubara frequencies shown [eV]
      siw       : (n_context,) raw Sigma_LL'(i wm), complex [eV]
      n_fit     : number of points actually used in the fit (resolved
                  from fit_wmax if given)
      fit_wm    : (n_curve,) dense frequency grid for the fit curve,
                  spanning [0, wm[n_context-1]] [eV]
      fit_re    : (n_curve,) fitted Re Sigma curve evaluated on fit_wm
      fit_im    : (n_curve,) fitted Im Sigma curve evaluated on fit_wm
      Sigma0    : fitted Re[Sigma(0+)]  [eV]
      Gamma     : fitted -Im[Sigma(0+)] [eV]
      Zinv      : fitted Z^-1 (only meaningful if L==Lp)
      diagonal  : bool, L==Lp
    """
    siw_work, _, N_total, _ = _get_siw_work(data, spin_average)
    wm = data.pos_iw

    if not (0 <= L < N_total and 0 <= Lp < N_total):
        raise ValueError("Orbital index out of range: (L,Lp)=({},{}), "
                          "N_total={}.".format(L, Lp, N_total))
    diagonal = (L == Lp)
    if data.offdiag:
        siw_ll = siw_work[L, Lp, :]
    else:
        if not diagonal:
            raise ValueError(
                "This is a diagonal-only run (qmc.offdiag=0); the "
                "off-diagonal element (L,Lp)=({},{}) was never read "
                "in (only siw, not siw-full, is available).".format(L, Lp))
        siw_ll = siw_work[L, :]

    n_fit_resolved = (_n_fit_from_wmax(wm, fit_wmax, beta=data.beta)
                      if fit_wmax is not None else n_fit)
    if n_fit_resolved > len(wm):
        raise ValueError("n_fit={} > Niw={}".format(n_fit_resolved, len(wm)))
    if n_fit_resolved <= fit_degree:
        raise ValueError(
            "n_fit={} <= fit_degree={}: a degree-{} polynomial needs at "
            "least {} points to be determined at all. Increase n_fit/"
            "fit_wmax or lower fit_degree.".format(
                n_fit_resolved, fit_degree, fit_degree, fit_degree + 1))
    if n_fit_resolved < 3 * (fit_degree + 1):
        logger.warning(
            "n_fit={} is less than 3x the number of fit parameters "
            "({} for a degree-{} polynomial); the fit may be poorly "
            "conditioned and offer little noise averaging.".format(
                n_fit_resolved, fit_degree + 1, fit_degree))

    Zinv_ll, Sigma0_ll, Gamma_ll = _fit_element(
        siw_ll, wm, n_fit_resolved, fit_degree, diagonal=diagonal)
    re_c, im_c = _polyfit_element(siw_ll, wm, n_fit_resolved, fit_degree)

    n_ctx = min(n_context, len(wm))
    fit_wm = np.linspace(0.0, wm[n_ctx - 1], n_curve)

    return dict(
        wm        = wm[:n_ctx],
        siw       = siw_ll[:n_ctx],
        n_fit     = n_fit_resolved,
        fit_wm    = fit_wm,
        fit_re    = np.polyval(re_c, fit_wm),
        fit_im    = np.polyval(im_c, fit_wm),
        Sigma0    = Sigma0_ll.real,
        Gamma     = Gamma_ll,
        Zinv      = Zinv_ll,
        diagonal  = diagonal,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Main extraction routine
# ─────────────────────────────────────────────────────────────────────────────

def extract_selfenergy(
    data: W2dynData,
    spin_average: Optional[bool] = None,
    n_fit: int = 10,
    fit_degree: int = 2,
    fit_wmax: Optional[float] = None,
    gamma_min: float = 1e-6,
) -> SelfEnergyParams:
    """
    Extract linearised QP parameters (Zinv_orb, Sigma0_orb, Gamma_orb) from
    a W2dynData object by polynomial fit to the lowest Matsubara frequencies.

    Parameters
    ----------
    data        : W2dynData from read_w2dyn().
    spin_average: True  = average spin channels (auto-True for 'para').
                  False = use spin-0 channel only (with a warning).
                  None  = auto-detect from data.magnetism.
    n_fit       : Number of lowest positive Matsubara points for the fit.
                  Ignored when fit_wmax is provided.
    fit_degree  : Polynomial degree, >= 1. Default 2 (quadratic): a linear
                  fit (1) is the minimal Fermi-liquid form but a quadratic
                  captures curvature within the fit window and is usually
                  preferable when statistics are good. Higher orders can
                  help further for runs with very good (e.g. high-statistics
                  QMC) data, at the cost of needing more points in the fit
                  window (n_fit) to stay well-conditioned -- a warning is
                  emitted if n_fit is less than ~3x the number of fit
                  parameters (fit_degree+1). Always inspect the fit
                  visually (structure/checkfit.py) before trusting a
                  higher-order extrapolation; see the accompanying tex,
                  sec:fitcheck, for why a higher-order *polynomial* (not a
                  spline) is still the right tool to reach for.
    fit_wmax    : If given [eV], use all Matsubara frequencies up to this
                  energy cutoff instead of a fixed n_fit count.  Useful for
                  comparing runs at different temperatures on a common energy
                  scale.  Takes priority over n_fit when set.
    gamma_min   : Floor for diagonal Gamma_LL [eV].  Applied with a warning
                  when the fit yields a negative (unphysical) value.

    Returns
    -------
    SelfEnergyParams
    """
    if spin_average is None:
        spin_average = (data.magnetism.lower() == 'para')
    logger.info("  spin_average={}  (magnetism='{}')".format(
        spin_average, data.magnetism))
    if not spin_average:
        logger.warning("spin_average=False: using spin-0 channel only.")

    wm  = data.pos_iw
    Niw = len(wm)

    # ── resolve n_fit from fit_wmax if provided ────────────────────────────
    if fit_wmax is not None:
        n_fit = _n_fit_from_wmax(wm, fit_wmax, beta=data.beta)
        logger.info("  fit_wmax={:.4f} eV  =>  n_fit={}".format(fit_wmax, n_fit))
    if n_fit > Niw:
        raise ValueError("n_fit={} > Niw={}".format(n_fit, Niw))
    if n_fit <= fit_degree:
        raise ValueError(
            "n_fit={} <= fit_degree={}: a degree-{} polynomial needs at "
            "least {} points to be determined at all, and considerably "
            "more to get any noise-averaging benefit out of the least-"
            "squares fit (that benefit is the whole point of preferring "
            "a polynomial fit over, e.g., a spline -- see the "
            "accompanying tex, sec:fitcheck). Increase n_fit/fit_wmax or "
            "lower fit_degree.".format(n_fit, fit_degree, fit_degree, fit_degree + 1))
    if n_fit < 3 * (fit_degree + 1):
        logger.warning(
            "n_fit={} is less than 3x the number of fit parameters "
            "({} for a degree-{} polynomial); the fit may be poorly "
            "conditioned and offer little noise averaging. Consider "
            "increasing n_fit/fit_wmax.".format(
                n_fit, fit_degree + 1, fit_degree))

    # ── spin averaging (shared with fit_diagnostics, see _get_siw_work) ────
    siw_work, dc_work, N_total, spin_average = _get_siw_work(data, spin_average)

    # ── allocate output matrices ───────────────────────────────────────────
    # Sigma0_orb : real-symmetric (N×N) [stored complex128, Im==0 always],
    #              holds Re[Sigma_LL'(0+)].
    # Gamma_orb  : real, symmetric, PSD (N×N), holds -Im[Sigma_LL'(0+)]
    #              for EVERY (L,L') pair (on and off the diagonal).
    # Zinv_orb   : real DIAGONAL (N×N), from the slope of Im[Sigma_LL].
    Zinv_orb   = np.eye(N_total, dtype=np.float64)
    Sigma0_orb = np.zeros((N_total, N_total), dtype=np.complex128)
    Gamma_orb  = np.zeros((N_total, N_total), dtype=np.float64)

    # ── fit ────────────────────────────────────────────────────────────────
    if data.offdiag:
        # Fit every (L, L') matrix element.
        # Sigma0_orb carries ONLY the real part for every pair (diagonal
        # and off-diagonal alike); it enters H_tilde and handles orbital
        # mixing without ever risking H_tilde's Hermiticity.
        # Gamma_orb is filled for every (L, L') pair; its off-diagonal
        # entries are still subject to Approximation 3 later (truncated
        # to the diagonal only AFTER rotation into the QP basis, see
        # wannier_dmft.py / tex sec:QPdiag) -- not here, at the orbital
        # level.
        for L in range(N_total):
            for Lp in range(N_total):
                Zinv_ll, S0_ll, G_ll = _fit_element(
                    siw_work[L, Lp, :], wm, n_fit, fit_degree, diagonal=(L==Lp))
                Sigma0_orb[L, Lp] = S0_ll
                Gamma_orb[L, Lp]  = G_ll
                if L == Lp:
                    Zinv_orb[L, L] = Zinv_ll

        # Symmetrise Sigma0_orb and Gamma_orb (suppresses fit noise between
        # the independently-fitted (L,L') and (L',L) entries).  Both are
        # real matrices (Sigma0_orb by construction, see _fit_element;
        # Gamma_orb because it is literally -Im of a single complex number
        # for each entry), so the physically correct symmetry relation is
        # plain transpose-symmetry, Sigma0_LL' ~= Sigma0_L'L and
        # Gamma_LL' ~= Gamma_L'L -- NOT Hermitian conjugation.  (Using
        # conjugate-transpose "Hermitisation" of the off-diagonal-complex
        # Sigma0_orb instead would be incorrect and would silently
        # discard the physical off-diagonal scattering information
        # whenever it is genuinely present -- see the module docstring.)
        max_asym_S0 = np.max(np.abs(Sigma0_orb - Sigma0_orb.T))
        if max_asym_S0 > 1e-3:
            logger.warning(
                "Sigma0_orb pre-symmetrisation asymmetry {:.2e} > 1e-3 eV; "
                "check siw-full convergence.".format(max_asym_S0))
        Sigma0_orb = 0.5 * (Sigma0_orb + Sigma0_orb.T)

        max_asym_G = np.max(np.abs(Gamma_orb - Gamma_orb.T))
        if max_asym_G > 1e-3:
            logger.warning(
                "Gamma_orb pre-symmetrisation asymmetry {:.2e} > 1e-3 eV; "
                "check siw-full convergence.".format(max_asym_G))
        Gamma_orb = 0.5 * (Gamma_orb + Gamma_orb.T)

        # Off-diagonal causality check: the full scattering-rate matrix
        # should be positive semi-definite (it is a Hermitian, and here
        # real-symmetric, broadening matrix; the diagonal Gamma_LL>=0
        # check below is the N=1 special case of this).
        evals_G = np.linalg.eigvalsh(Gamma_orb)
        if np.any(evals_G < -1e-6):
            logger.warning(
                "Gamma_orb has negative eigenvalues (min={:.3e} eV); "
                "the off-diagonal scattering matrix is not positive "
                "semi-definite.  Check siw-full convergence / "
                "N_fit.".format(evals_G.min()))

    else:
        # Diagonal run: only N diagonal elements; all off-diagonal are 0.
        for L in range(N_total):
            Zinv_ll, S0_ll, G_ll = _fit_element(
                siw_work[L, :], wm, n_fit, fit_degree, diagonal=True)
            Zinv_orb[L, L]   = Zinv_ll
            Sigma0_orb[L, L] = S0_ll
            Gamma_orb[L, L]  = G_ll

    # ── sanity checks ──────────────────────────────────────────────────────
    diag_Zinv  = np.diag(Zinv_orb)
    diag_Z     = 1.0 / diag_Zinv     # Z itself: more physically intuitive than Z^-1
    diag_Gamma = np.diag(Gamma_orb)
    diag_Sigma0 = np.diag(Sigma0_orb).real

    if np.any(diag_Zinv < 1.0 - 1e-6):
        logger.warning("Zinv diagonal has entries < 1: {}.  "
                       "Check DMFT convergence.".format(diag_Zinv))

    # Display layer only (does not affect any internal indexing): orbitals
    # are numbered 1..N here and in fit_diagnostics/checkfit.py plots, the
    # usual physics convention, even though everything internally (arrays,
    # the --elements CLI argument of checkfit.py, etc.) stays 0-indexed.
    logger.info("  per-orbital fit results (orbital index is 1-based):")
    for L in range(N_total):
        logger.info(
            "    orbital {:2d}:  Z={:.5f}  (Z^-1={:.5f})   "
            "Sigma0={:+.5f} eV   Gamma={:.5f} eV".format(
                L + 1, diag_Z[L], diag_Zinv[L], diag_Sigma0[L], diag_Gamma[L]))

    n_neg = np.sum(diag_Gamma < 0)
    if n_neg > 0:
        logger.warning("{} diagonal Gamma < 0; applying floor {:.1e} eV.  "
                       "Check convergence.".format(n_neg, gamma_min))
        for L in range(N_total):
            if Gamma_orb[L, L] < gamma_min:
                Gamma_orb[L, L] = gamma_min

    evals = np.linalg.eigvalsh(Zinv_orb)
    if np.any(evals <= 0):
        raise ValueError(
            "Zinv_orb not positive-definite (min eigenvalue {:.6f}).  "
            "Cannot proceed with GEP.".format(evals.min()))

    return SelfEnergyParams(
        temperature = data.temperature,
        beta        = data.beta,
        mu          = data.mu,
        n_orb       = N_total,
        spin_average= spin_average,
        Zinv_orb    = Zinv_orb,
        Sigma0_orb  = Sigma0_orb,
        Gamma_orb   = Gamma_orb,
        dc_orb      = dc_work,
        n_fit       = n_fit,
        fit_degree  = fit_degree,
        fit_wmax    = fit_wmax,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Multi-file convenience
# ─────────────────────────────────────────────────────────────────────────────

def read_and_sort_temperatures(
    paths: list,
    spin_average: Optional[bool] = None,
    n_fit: int = 10,
    fit_degree: int = 2,
    fit_wmax: Optional[float] = None,
    gamma_min: float = 1e-6,
) -> list:
    """
    Read multiple w2dynamics HDF5 files, extract QP parameters, return
    sorted by ascending temperature.

    When fit_wmax is provided it takes priority over n_fit for every file,
    ensuring the fit uses the same energy window across temperatures
    (though different numbers of Matsubara points due to different beta).
    """
    if not paths:
        raise ValueError("No HDF5 paths provided.")

    params_list = []
    for p in paths:
        data   = read_w2dyn(p)
        params = extract_selfenergy(data,
                                    spin_average=spin_average,
                                    n_fit=n_fit,
                                    fit_degree=fit_degree,
                                    fit_wmax=fit_wmax,
                                    gamma_min=gamma_min)
        params_list.append(params)
        logger.info("  T={:.2f} K  n_fit={}  from {}".format(
            params.temperature, params.n_fit, p))

    params_list.sort(key=lambda p: p.temperature)
    logger.info("Temperature series (K): {}".format(
        ["{:.2f}".format(p.temperature) for p in params_list]))
    return params_list
