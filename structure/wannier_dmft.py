#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
structure/wannier_dmft.py
==========================
Provides the `computeDMFT` method to be grafted onto `Wannier90Calculation`.

Kept separate from wannier.py to minimise the diff; bound at import time in
lwann:

    from structure.wannier_dmft import computeDMFT as _computeDMFT
    Wannier90Calculation.computeDMFT = _computeDMFT

Algorithm (see dmft_lintretrace_interface.tex)
-------------------------------------------------
For each k-point on the reducible BZ mesh:

  Step 2  FT in orbital basis: H(k), hv(k) [+Peierls], hc(k)
  Step 3  Form H_tilde(k) = H(k) + Sigma0_orb + diag(dc) - mu*I in orbital
          basis (Sigma0_orb is real-symmetric for every orbital pair, on
          and off the diagonal -- see structure/w2dyn.py docstring -- so
          H_tilde stays exactly Hermitian).
          Solve scipy.linalg.eigh(H_tilde, Zinv_orb) -> E_QP(k), W(k).
          Also diagonalise bare H(k) -> eps_bare(k) for energy-file ordering.
  Step 4  Rotate hv and hc with W in a single pass:
              vk_QP = W† hv(k) W,   ck_QP = W† hc(k) W
          Compute optical moments exactly as in the lwann reducible-grid path.
  Step 5  Extract per-band QP parameters (see _per_k_dmft):
              gamma_QP[n] = [W† Gamma_orb W]_nn   (real; Gamma_orb may now
                            carry off-diagonal orbital entries too, see
                            structure/w2dyn.py)
              Z_QP[n]     = [W† W]_nn             (real; NOT
                            1/[W† Zinv_orb W]_nn, see tex sec:QPdiag)
              bandshift[n] = E_QP[n]/Z_QP[n] - eps_bare[n] + mu_dmft
                            (NOT simply E_QP[n] - eps_bare[n]; the /Z_QP
                            and the restoration of mu_dmft are both
                            required because LinReTraCe's own read-in
                            formula band = band_file + band_shift is
                            applied BEFORE the Z-multiplication, see tex
                            eq:lrtc_enrgy/eq:bandshift_final).
  Step 6  Before writing to the scattering HDF5: gamma_write[n] =
          gamma_QP[n] / Z_QP[n].  This division is REQUIRED because
          LinReTraCe's Fortran core (src_linretrace/input.f90)
          performs its OWN multiplication `sct%gam = sct%gam * sct%zqp`
          on read-in.  gamma_QP already equals the fully Z-renormalised
          physical rate (tex sec:QPdiag); writing it unmodified to
          'scatrate' would let LinReTraCe multiply it by Z_QP a SECOND
          time, silently double-suppressing the scattering rate (e.g.
          by Z^2 instead of Z in the single-orbital case).  Z_QP and
          bandshift are written as computed in Step 5, unmodified.
  Step 7  Off-diagonal scattering diagnostic (Approximation 3, tex
          sec:validity): at every k, _per_k_dmft also returns an
          OffdiagDiag bundling the Frobenius norms of the diagonal and
          off-diagonal blocks of W† Gamma_orb W (the latter is what
          Approximation 3 discards) -- both raw (unweighted) and
          weighted by transport relevance: a Fermi-window factor
          (suppresses pairs of QP states far from mu, in units of k_B T;
          the same formula handles insulators automatically, since f->0/1
          flattens deep in either band) times a Lorentzian coherence
          factor (suppresses pairs whose energy splitting exceeds their
          combined linewidth, i.e. resolved/incoherent peaks).  After
          each (T, spin) k-loop, _report_offdiag_scattering_diagnostic
          logs BZ-aggregate ratios (raw and weighted) plus the worst
          individual k-points by the weighted ratio (location, dominant
          QP band pair, and the energies/linewidths that explain the
          weighting), so the validity of Approximation 3 for ORDINARY
          (non-anomalous) DC transport can be assessed without a
          separate post-processing pass.

Conventions (matching lwann/computeHamiltonian exactly)
-------------------------------------------------------
  - hv(k) shape  : (norb, norb, 3)   index order (L, L', alpha)
  - hc(k) shape  : (norb, norb, 6)   index order (L, L', dir)
    where dir enumerates (xx,yy,zz,xy,xz,yz)
  - W shape       : (norb, nbands)    W[:,n] is the n-th QP eigenvector
  - Peierls correction (matching lwann exactly):
      distances   = plist[:,None,:] - plist[None,:,:]  (tau_L - tau_Lp, fractional)
      ri_minus_rj = einsum('id,abi->abd', rvec, distances)  (Cartesian)
      correction  = -1j * H(k)[:,:,None] * ri_minus_rj
  - Non-ortho optical elements shape (9):
      indices 0..5 : real parts of (v*v)  for pairs (xx,yy,zz,xy,xz,yz)
      indices 6..8 : imaginary parts of (v*v) for diagonal pairs (xx,yy,zz)
      (this matches lwann exactly)

Memory: O(norb^2) working arrays per k-point.
Threading: BLAS/OpenMP from OMP_NUM_THREADS is inherited automatically.
"""

from __future__ import annotations

import logging
import os
from typing import List, NamedTuple, Tuple

import numpy as np
import scipy.linalg

from structure.w2dyn import SelfEnergyParams
from scattering.fullscattering import FullScattering

logger = logging.getLogger(__name__)

_KB_EV = 8.61733034e-5   # Boltzmann constant [eV / K]


# ─────────────────────────────────────────────────────────────────────────────
# Levi-Civita symbol (computed once at module load)
# ─────────────────────────────────────────────────────────────────────────────

def _levi_civita():
    eps = np.zeros((3, 3, 3), dtype=np.float64)
    for i, j, k in [(0,1,2),(1,2,0),(2,0,1)]:
        eps[i,j,k] =  1.0
        eps[k,j,i] = -1.0
    return eps

_LEV = _levi_civita()


# ─────────────────────────────────────────────────────────────────────────────
# Off-diagonal scattering diagnostic bundle (Approximation 3, tex sec:validity)
# ─────────────────────────────────────────────────────────────────────────────

class OffdiagDiag(NamedTuple):
    """
    Per-k bundle quantifying how much off-diagonal QP-basis scattering
    Approximation 3 discards, in both a raw (unweighted) and a transport-
    relevance-weighted sense.  See _per_k_dmft and tex sec:validity item 5.

    Raw (unweighted) -- purely about the mathematical size of the
    discarded matrix, e.g. useful for spotting fit/convergence issues
    independent of whether the relevant states matter for transport:
      diag_norm, offdiag_norm : float   ||.||_F of diag/off-diag blocks
      offdiag_max             : float   largest |off-diagonal entry| [eV]
      offdiag_pair            : (int,int)  0-based QP band indices of that entry

    Transport-relevance-weighted -- the same quantities after multiplying
    every (n,m) matrix element by a weight in (0,1] that captures two
    physically-motivated reasons a large off-diagonal entry might still
    not matter for ORDINARY (non-anomalous) DC transport (see _per_k_dmft
    for the formula): both QP states being close to mu (in units of k_B T;
    the SAME formula also handles insulators, since it automatically
    suppresses states deep in either band), and their energy splitting
    not exceeding their combined linewidth (Lorentzian coherence factor).
      diag_norm_w, offdiag_norm_w : float  weighted ||.||_F
      offdiag_max_w               : float  largest weighted |entry| [eV]
      offdiag_pair_w              : (int,int)  0-based QP band indices
      offdiag_pair_w_E            : (float,float)  E_n-mu, E_m-mu of that
                                     pair [eV], for interpreting the weight
      offdiag_pair_w_gamma        : (float,float)  gamma_n, gamma_m of
                                     that pair [eV], ditto
    """
    diag_norm: float
    offdiag_norm: float
    offdiag_max: float
    offdiag_pair: Tuple[int, int]
    diag_norm_w: float
    offdiag_norm_w: float
    offdiag_max_w: float
    offdiag_pair_w: Tuple[int, int]
    offdiag_pair_w_E: Tuple[float, float]
    offdiag_pair_w_gamma: Tuple[float, float]


# ─────────────────────────────────────────────────────────────────────────────
# Per-k DMFT pipeline
# ─────────────────────────────────────────────────────────────────────────────

def _per_k_dmft(hk_orb, hvk_orb, hck_orb,
                Zinv_orb, Sigma0dc_orb, Gamma_orb,
                mu_dmft, beta_dmft, ortho):
    """
    Full DMFT renormalisation for a single k-point.

    Parameters
    ----------
    hk_orb        : (norb, norb) complex  –  bare H(k) in orbital basis
    hvk_orb       : (norb, norb, 3) complex  –  velocity H (Peierls applied)
    hck_orb       : (norb, norb, 6) complex  –  curvature H
    Zinv_orb      : (norb, norb) float    –  inverse QP weight (k-independent)
    Sigma0dc_orb  : (norb, norb) complex  –  Sigma0_orb + diag(dc_orb), full
                     real-symmetric matrix (Im==0 by construction, see
                     structure/w2dyn.py), diagonal-only for offdiag=False
                     runs but carries off-diagonal orbital mixing for
                     offdiag=True
    Gamma_orb     : (norb, norb) float    –  orbital scattering-rate matrix
                     (k-independent).  Diagonal-only for offdiag=False
                     runs; for offdiag=True runs may carry off-diagonal
                     orbital entries too (symmetric, PSD, see
                     structure/w2dyn.py).  Only its diagonal in the QP
                     basis survives projection below (Approximation 3).
    mu_dmft       : float                 –  DMFT chemical potential [eV]
    beta_dmft     : float                 –  DMFT inverse temperature [eV^-1],
                     used only for the transport-relevance weighting of
                     the off-diagonal scattering diagnostic below (it does
                     NOT otherwise enter the GEP/QP extraction at all).
    ortho         : bool                  –  orthogonal unit cell flag

    Returns
    -------
    eps_bare  : (norb,) float   bare KS eigenvalues (ascending)
    E_QP      : (norb,) float   QP energies (ascending)
    gamma_QP  : (norb,) float   per-band scattering rate [eV], already
                fully Z-renormalised (tex sec:QPdiag); NOT what gets
                written to 'scatrate' directly -- see computeDMFT, which
                divides by Z_QP first to compensate LinReTraCe's own
                read-in convention.
    Z_QP      : (norb,) float   per-band QP weight, [W^dagger W]_nn
    bandshift : (norb,) float   E_QP/Z_QP - eps_bare + mu_dmft [eV]
                (NOT simply E_QP - eps_bare; see module docstring Step 5)
    vel2diag  : (norb, 3 or 9) float   diagonal optical moments
    Bvel2diag : (norb, 3, 3, 3) complex  diagonal B-field optical moments
    vel2      : (norb, norb, 3 or 9) float  full optical moments
    Bvel2     : (norb, norb, 3, 3, 3) complex  full B-field moments
    diag      : OffdiagDiag   off-diagonal scattering diagnostic bundle,
                see its docstring and tex sec:validity item 5.
    """
    norb = hk_orb.shape[0]

    # ── bare KS eigenvalues (energy-file ordering reference) ─────────────
    eps_bare, _ = np.linalg.eigh(hk_orb)
    eps_bare = np.sort(eps_bare.real)

    # ── effective Hermitian matrix in orbital basis (Step 3) ─────────────
    # Full matrix addition: includes orbital off-diagonal self-energy
    # (qmc.offdiag=1 / siw-full runs) automatically.  For the diagonal
    # case Sigma0dc_orb is diagonal, so this reduces exactly to the
    # earlier diagonal-only update.
    diag_idx = np.arange(norb)
    H_tilde = hk_orb + Sigma0dc_orb
    H_tilde[diag_idx, diag_idx] -= mu_dmft

    # ── generalised Hermitian eigenvalue problem (or 1-band fast path) ───
    if norb == 1:
        # 1-band scalar case: bypass scipy.linalg.eigh completely.
        # H_tilde is 1×1 complex: H = [[h]], Zinv = [[z]], z > 0.
        # Generalised eigenproblem: h*w = E*z*w => E = h/z = Z*h.
        # W = [[1/sqrt(z)]] from the normalisation w†*z*w = 1.
        h = H_tilde[0, 0].real   # Hermitian 1×1 diagonal is exactly real
        z = Zinv_orb[0, 0]
        Z_val = 1.0 / z
        E_QP   = np.array([h * Z_val])
        W      = np.array([[1.0 / np.sqrt(z)]])
        Wdag   = W.conj().T
    else:
        # General N-band case: scipy.linalg.eigh(A, B) solves A v = λ B v.
        E_QP, W = scipy.linalg.eigh(H_tilde, Zinv_orb)
        E_QP = E_QP.real     # already ascending from eigh
        Wdag = W.conj().T

    WdagGW   = Wdag @ Gamma_orb @ W
    gamma_QP = np.diag(WdagGW).real

    # Physical QP weight: residue of G(k,ω) at pole E_n.
    # From the generalised orthonormality W† Zinv W = I the eigh eigenvectors
    # are normalised as ||W[:,n]||²_{Zinv} = 1, NOT ||W[:,n]||² = 1.
    # The actual spectral weight is the squared norm in the standard inner product:
    #   Z̃_n = (W† W)_{nn} = sum_L |W_{Ln}|²
    # This satisfies the sum rule: sum_n Z̃_n = Tr(W† W) = Tr(Z^{orb}).
    # (Computed here, ahead of where it was previously needed, because the
    # transport-relevance weighting below also needs E_QP -- already have
    # it -- and nothing else new; Z_QP itself doesn't enter the weighting.)
    Z_QP = np.diag(Wdag @ W).real

    # ── off-diagonal scattering diagnostic (Approximation 3) ─────────────
    # WdagGW is the FULL scattering-rate matrix already rotated into the
    # QP/band basis; gamma_QP above keeps only its diagonal (Approximation
    # 3, tex sec:QPdiag).  Quantify, at this k-point, how much is being
    # discarded: the Frobenius norm of the off-diagonal block relative to
    # the diagonal block.  Essentially free -- WdagGW is already computed.
    diag_norm = np.linalg.norm(gamma_QP)               # ||Gamma~^d(k)||_F
    offdiag   = WdagGW.copy()
    offdiag[diag_idx, diag_idx] = 0.0                   # zero the diagonal
    offdiag_norm = np.linalg.norm(offdiag)              # ||Gamma~^od(k)||_F
    abs_offdiag  = np.abs(offdiag)
    if norb > 1:
        n_dom, m_dom = np.unravel_index(np.argmax(abs_offdiag), abs_offdiag.shape)
        offdiag_max  = abs_offdiag[n_dom, m_dom]
    else:
        n_dom, m_dom, offdiag_max = 0, 0, 0.0           # no off-diagonal possible

    # ── transport-relevance weighting of the same diagnostic ──────────────
    # A large |Gamma~^od_nm(k)| does not automatically matter for ORDINARY
    # (non-anomalous) DC transport.  Two physically-motivated suppression
    # mechanisms, combined into one dimensionless weight w_nm in (0,1]:
    #
    #  (i) Fermi window.  The Boltzmann/Kubo DC conductivity kernel is
    #      weighted by -df/dE, sharply peaked at E=mu with width ~k_B T.
    #      Per-band weight, normalised to 1 exactly at E_n=mu:
    #          fermi_w(E_n) = 1/cosh^2( beta*(E_n-mu)/2 ) = 4 f(1-f)
    #      This SAME formula, with no separate branch, also correctly
    #      handles insulators: deep in the valence (f->1) or conduction
    #      (f->0) band fermi_w->0, so only near-band-edge states (those
    #      within a few k_B T of mu, wherever mu sits) keep any weight --
    #      exactly the qualitative behaviour wanted, automatically, with
    #      no metal/insulator switch.  (In a wide-gap insulator at low T
    #      this can suppress the weighted ratio towards 0 everywhere, even
    #      at the band edges; that is not a flaw, it is the correct
    #      statement that essentially no DC transport -- and hence no
    #      transport correction from neglected interband scattering --
    #      is happening at all in that regime.)
    #  (ii) Coherence.  Two QP "peaks" of width gamma_n, gamma_m separated
    #      by |E_n-E_m| only meaningfully overlap/mix for transport if
    #      that splitting does not exceed their combined width; modelled
    #      with a standard Lorentzian overlap factor,
    #          coh_w(n,m) = (gamma_n+gamma_m)^2 / [(E_n-E_m)^2+(gamma_n+gamma_m)^2]
    #      which is exactly 1 for degenerate/strongly-broadened pairs and
    #      falls off like [(gamma_n+gamma_m)/(E_n-E_m)]^2 for well-resolved
    #      ones.  This is a qualitative, order-of-magnitude weighting (the
    #      literal vertex-correction prefactor would require redoing the
    #      Kubo bubble in this basis), not a quantitatively exact one --
    #      see tex sec:validity for the caveat.
    # w_nm = sqrt(fermi_w(E_n)*fermi_w(E_m)) * coh_w(n,m); reduces to
    # fermi_w(E_n) exactly on the diagonal (coh_w=1 there), so the weighted
    # diagonal norm below uses the same single-band Fermi-window weight one
    # would naively expect.
    fermi_w = 1.0 / np.cosh(np.clip(beta_dmft * (E_QP - mu_dmft) / 2.0, -300, 300))**2
    gam_sum = gamma_QP[:, None] + gamma_QP[None, :]
    dE      = E_QP[:, None] - E_QP[None, :]
    coh_w   = gam_sum**2 / (dE**2 + gam_sum**2 + 1e-300)
    pair_w  = np.sqrt(fermi_w[:, None] * fermi_w[None, :]) * coh_w

    WdagGW_w = pair_w * WdagGW
    diag_norm_w = np.linalg.norm(np.diag(WdagGW_w).real)
    offdiag_w   = WdagGW_w.copy()
    offdiag_w[diag_idx, diag_idx] = 0.0
    offdiag_norm_w = np.linalg.norm(offdiag_w)
    abs_offdiag_w  = np.abs(offdiag_w)
    if norb > 1:
        n_dom_w, m_dom_w = np.unravel_index(np.argmax(abs_offdiag_w), abs_offdiag_w.shape)
        offdiag_max_w = abs_offdiag_w[n_dom_w, m_dom_w]
    else:
        n_dom_w, m_dom_w, offdiag_max_w = 0, 0, 0.0

    diag = OffdiagDiag(
        diag_norm=diag_norm, offdiag_norm=offdiag_norm,
        offdiag_max=offdiag_max, offdiag_pair=(int(n_dom), int(m_dom)),
        diag_norm_w=diag_norm_w, offdiag_norm_w=offdiag_norm_w,
        offdiag_max_w=offdiag_max_w, offdiag_pair_w=(int(n_dom_w), int(m_dom_w)),
        offdiag_pair_w_E=(float(E_QP[n_dom_w] - mu_dmft),
                          float(E_QP[m_dom_w] - mu_dmft)),
        offdiag_pair_w_gamma=(float(gamma_QP[n_dom_w]), float(gamma_QP[m_dom_w])),
    )

    # ── band shift: must invert LinReTraCe's OWN reconstruction formula ──
    # LinReTraCe core (src_linretrace/input.f90) does, at read-in:
    #     edisp%band  = edisp%band_file + edisp%band_shift          (Eq. A)
    # and then, throughout response.F90/root.F90:
    #     enrgy(n,k)  = sct%zqp(n,k) * (edisp%band(n,k) - mu)        (Eq. B)
    # where edisp%band_file == eps_bare (the energy-file content) and mu is
    # whatever chemical potential LinReTraCe is evaluating at that moment.
    #
    # We want enrgy(n,k) to reproduce the true DMFT QP energy measured from
    # mu_dmft, i.e. enrgy(n,k) == E_QP[n] when mu == mu_dmft.  Substituting
    # Eq. A into Eq. B and solving for bandshift:
    #
    #     Z_QP[n] * (eps_bare[n] + bandshift[n] - mu_dmft) = E_QP[n]
    #     bandshift[n] = E_QP[n] / Z_QP[n] - eps_bare[n] + mu_dmft
    #
    # NOTE: this is NOT simply (E_QP - eps_bare): omitting the 1/Z_QP
    # factor or the mu_dmft restoration would silently corrupt both the
    # band positions and (through the implicit T,mu dependence of Z) the
    # apparent bandwidth at every temperature.  The two formulas only
    # coincide in the special case Z_QP=1 (no Z-renormalisation at all),
    # so this is easy to miss in a 1-orbital-only test; see
    # dmft_lintretrace_interface.tex, sec:twofile, for the full derivation.
    bandshift = E_QP / Z_QP - eps_bare + mu_dmft

    # ── rotate hv and hc by W in a single pass (Step 4) ──────────────────
    # vk_QP[n, m, alpha] = sum_{L,L'} W*[L,n] hvk[L,L',alpha] W[L',m]
    #                     = (Wdag @ hvk[:,:,a] @ W)[n,m]
    # Using einsum over the orbital indices, keeping alpha as a free index:
    vk_QP = np.einsum('ni,ija,jm->nma', Wdag, hvk_orb, W)   # (norb,norb,3)

    # ck_QP[n, m, dir] = (Wdag @ hck[:,:,d] @ W)[n,m]
    ck_QP = np.einsum('ni,ija,jm->nma', Wdag, hck_orb, W)   # (norb,norb,6)

    # ── optical elements (matching lwann reducible-grid path exactly) ─────
    vk_QP_conj = np.conjugate(vk_QP)

    # expand curvature to full 3×3 symmetric matrix form
    curmat = np.zeros((norb, norb, 3, 3), dtype=np.complex128)
    curmat[:, :, [0,1,2,0,0,1], [0,1,2,1,2,2]] = ck_QP
    curmat[:, :, [1,2,2], [0,0,1]] = curmat[:, :, [0,0,1], [1,2,2]]

    # velocity-squared optical elements: v*_a * v_b for 6 pair combinations
    vel2_raw = vk_QP_conj[:, :, [0,1,2,0,0,1]] * vk_QP[:, :, [0,1,2,1,2,2]]
    # shape: (norb, norb, 6)

    if ortho:
        vel2 = vel2_raw[:, :, :3].real          # (norb, norb, 3)
    else:
        temp = vel2_raw.copy()
        vel2 = np.empty((norb, norb, 9), dtype=np.float64)
        vel2[:, :, :6] = temp.real
        vel2[:, :, 6:] = temp[:, :, :3].imag   # imag parts of diag pairs (xx,yy,zz)

    vel2diag  = vel2[diag_idx, diag_idx, :]     # (norb, 3 or 9)

    # B-field optical moments: epsilon_cij * v*_a * v_i * c_bj
    Bvel2     = np.einsum('cij,nma,nmi,nmbj->nmabc',
                          _LEV, vk_QP_conj, vk_QP, curmat)  # (norb,norb,3,3,3)
    Bvel2diag = Bvel2[diag_idx, diag_idx, :, :, :]           # (norb, 3, 3, 3)

    return (eps_bare, E_QP, gamma_QP, Z_QP, bandshift,
            vel2diag, Bvel2diag, vel2, Bvel2, diag)


# ─────────────────────────────────────────────────────────────────────────────
# Off-diagonal scattering diagnostic (Approximation 3, tex sec:validity)
# ─────────────────────────────────────────────────────────────────────────────

def _report_offdiag_scattering_diagnostic(
        diags: List["OffdiagDiag"], kpoints, temperature, ispin, nspins, n_worst=5):
    """
    Log a summary of how much off-diagonal QP-basis scattering
    Approximation 3 discards, for one (temperature, spin) combination,
    both raw (unweighted) and weighted by transport relevance.

    Three numbers per version (raw / weighted), as discussed in tex
    sec:validity:

      1. BZ-aggregate ratio, RMS-combined over k-points:
             R_agg = sqrt(sum_k ||Gamma~^od(k)||_F^2)
                     / sqrt(sum_k ||Gamma~^d(k)||_F^2)
         (using the weighted diag/off-diag norms for the weighted
         version of R_agg, so numerator and denominator are weighted
         consistently). This is dominated by k-points with large
         absolute scattering, so a few small-denominator outliers
         cannot blow it up artificially -- unlike the plain mean of
         the per-k ratio, also reported for comparison; a large gap
         between the two flags exactly this kind of outlier.
      2. A localised metric (weighted version only, since it is the
         actionable one): the n_worst individual k-points with the
         largest weighted per-k ratio, together with the dominant
         (largest-magnitude, weighted) off-diagonal QP band pair, its
         distance from mu, and its combined linewidth -- i.e. not just
         "how bad" but "why" (close to mu? broad? resolved or not?).

    Parameters
    ----------
    diags : list[OffdiagDiag], length nkp (one per k-point, see
            _per_k_dmft)
    kpoints : (nkp, 3) float   fractional k-point coordinates
    temperature : float [K]
    ispin, nspins : int   (ispin is 0-based; reported 1-based, like the
                  existing "spin {}/{}" log line elsewhere in this module)
    n_worst : int   how many worst k-points to list (default 5)
    """
    nkp = len(diags)
    eps = 1e-300   # guards against an all-zero Gamma_orb (no scattering at all)

    diag_norm_k        = np.array([d.diag_norm for d in diags])
    offdiag_norm_k      = np.array([d.offdiag_norm for d in diags])
    diag_norm_w_k       = np.array([d.diag_norm_w for d in diags])
    offdiag_norm_w_k    = np.array([d.offdiag_norm_w for d in diags])
    offdiag_max_w_k     = np.array([d.offdiag_max_w for d in diags])
    offdiag_pair_w_k    = np.array([d.offdiag_pair_w for d in diags])
    offdiag_pair_w_E_k  = np.array([d.offdiag_pair_w_E for d in diags])
    offdiag_pair_w_gam_k= np.array([d.offdiag_pair_w_gamma for d in diags])

    R_agg   = np.sqrt(np.sum(offdiag_norm_k**2))   / max(np.sqrt(np.sum(diag_norm_k**2)), eps)
    R_mean  = float(np.mean(offdiag_norm_k   / np.maximum(diag_norm_k,   eps)))
    R_agg_w = np.sqrt(np.sum(offdiag_norm_w_k**2)) / max(np.sqrt(np.sum(diag_norm_w_k**2)), eps)
    per_k_ratio_w = offdiag_norm_w_k / np.maximum(diag_norm_w_k, eps)
    R_mean_w = float(np.mean(per_k_ratio_w))

    logger.info(
        "  off-diagonal scattering diagnostic (Approximation 3), "
        "T={:.2f} K, spin {}/{}:".format(temperature, ispin + 1, nspins))
    logger.info(
        "    raw       BZ-aggregate ratio = {:.4f}   mean per-k ratio = {:.4f}  "
        "(RMS over {} k-points; purely the size of what is discarded, "
        "ignoring transport relevance)".format(R_agg, R_mean, nkp))
    logger.info(
        "    weighted  BZ-aggregate ratio = {:.4f}   mean per-k ratio = {:.4f}  "
        "(Fermi-window x coherence weighted -- this is the one that "
        "matters for DC transport)".format(R_agg_w, R_mean_w))

    if nkp > 1 and np.any(offdiag_norm_w_k > 0):
        n_show = min(n_worst, nkp)
        worst  = np.argsort(per_k_ratio_w)[::-1][:n_show]
        logger.info("    worst {} k-points by WEIGHTED per-k ratio "
                    "(band indices below are 1-based):".format(n_show))
        for rank, ik in enumerate(worst):
            kx, ky, kz = kpoints[ik]
            n_dom, m_dom = offdiag_pair_w_k[ik]
            En, Em = offdiag_pair_w_E_k[ik]
            gn, gm = offdiag_pair_w_gam_k[ik]
            logger.info(
                "      #{:<2d} ratio={:.4f}  k=({:+.4f},{:+.4f},{:+.4f})  "
                "(ik={})   dominant pair: bands ({},{})  "
                "|Gamma~^od_w|={:.5f} eV   "
                "[E-mu=({:+.4f},{:+.4f}) eV, gamma=({:.4f},{:.4f}) eV]".format(
                    rank + 1, per_k_ratio_w[ik], kx, ky, kz, ik + 1,
                    n_dom + 1, m_dom + 1, offdiag_max_w_k[ik], En, Em, gn, gm))


# ─────────────────────────────────────────────────────────────────────────────
# Main method  (bound onto Wannier90Calculation in lwann)
# ─────────────────────────────────────────────────────────────────────────────

def computeDMFT(self,
                se_params_list: List[SelfEnergyParams],
                energy_hdf5: str,
                output: str = 'lrtc-dmft-scat.hdf5',
                peierlscorrection: bool = True,
                n_wann_check: bool = True):
    """
    Generate the LinReTraCe scattering HDF5 from a list of SelfEnergyParams.

    Must be called after readData() (and optionally expandKmesh()).
    Does NOT require computeHamiltonian() to have been called; the energy HDF5
    written by a preceding plain lwann run is read to initialise FullScattering.

    Parameters
    ----------
    se_params_list : list of SelfEnergyParams
        One entry per temperature, sorted ascending (from read_and_sort_temperatures).
    energy_hdf5 : str
        Path to the existing LinReTraCe energy HDF5.
    output : str
        Path for the output scattering HDF5.
    peierlscorrection : bool
        Apply inter-orbital Peierls correction to velocity. Default True.
    n_wann_check : bool
        Assert nproj == n_orb from DMFT. Default True.
    """
    logger.info("="*60)
    logger.info("computeDMFT: generating DMFT scattering file")
    logger.info("  energy file : {}".format(energy_hdf5))
    logger.info("  output      : {}".format(output))
    logger.info("  k-mesh      : {}x{}x{} = {} k-points".format(
        self.nkx, self.nky, self.nkz, self.nkp))

    if not se_params_list:
        raise ValueError("se_params_list is empty.")
    nt = len(se_params_list)

    # ── orbital count consistency check ───────────────────────────────────
    n_orb_wann = self.nproj
    n_orb_dmft = se_params_list[0].n_orb
    if n_wann_check and n_orb_wann != n_orb_dmft:
        raise ValueError(
            "Wannier90 nproj={} != DMFT n_orb={}. "
            "The self-energy must cover the full Wannier basis.".format(
                n_orb_wann, n_orb_dmft))

    # ── spin-average check ────────────────────────────────────────────────
    for sp in se_params_list:
        if not sp.spin_average:
            raise NotImplementedError(
                "computeDMFT currently requires spin_average=True. "
                "Full spin-resolved output is not yet implemented.")

    # ── initialise FullScattering from the energy file ────────────────────
    scatobj = FullScattering(energy_hdf5)
    spins, nkp, nbands = scatobj.getDependencies()
    if nkp != self.nkp:
        raise ValueError(
            "k-mesh in energy HDF5 ({} k-points) differs from current "
            "Wannier90 k-mesh ({} k-points).".format(nkp, self.nkp))
    if nbands != n_orb_wann:
        raise ValueError(
            "Band count in energy HDF5 ({}) != nproj ({}).".format(
                nbands, n_orb_wann))
    if spins != self.spins:
        raise ValueError(
            "Spin count in energy HDF5 ({}) != Wannier90Calculation spin "
            "count ({}). The energy file (.bands/ispin) and this "
            "Wannier90 run must come from the same hr.dat / calctype; "
            "gamma_out/Z_out/bandshift_out below are shaped by the "
            "former but the k-loop below iterates over the latter, "
            "which would otherwise silently mismatch.".format(
                spins, self.spins))

    # ── temperature axis ──────────────────────────────────────────────────
    temps = np.array([sp.temperature for sp in se_params_list], dtype=np.float64)
    tmin, tmax = float(temps[0]), float(temps[-1])
    # defineTemperatures creates a linspace; we override with exact DMFT temps
    scatobj.defineTemperatures(tmin=tmin, tmax=tmax, nt=nt, tlog=False)
    scatobj.temps = temps
    scatobj.betas = 1.0 / (temps * _KB_EV)
    scatobj.mus   = np.full(nt, scatobj.mudft, dtype=np.float64)
    scatobj.tmin  = tmin
    scatobj.tmax  = tmax

    # ── pre-allocate output arrays ────────────────────────────────────────
    gamma_out     = np.zeros((nt, spins, nkp, nbands), dtype=np.float64)
    Z_out         = np.ones ((nt, spins, nkp, nbands), dtype=np.float64)
    bandshift_out = np.zeros((nt, spins, nkp, nbands), dtype=np.float64)

    # ── k-independent prefactors for FT ──────────────────────────────────
    # prefactor_r[d, r]  = R_Cartesian[r, d]  shape (3, nrpts)
    prefactor_r  = np.einsum('id,ri->dr', self.rvec, self.rlist)

    # prefactor_r2[dir, r] = R_d1 * R_d2  for 6 direction pairs (xx,yy,zz,xy,xz,yz)
    prefactor_r2 = np.zeros((6, self.nrp), dtype=np.float64)
    for idir, i, j in zip(range(6), [0,1,2,0,0,1], [0,1,2,1,2,2]):
        prefactor_r2[idir, :] = prefactor_r[i, :] * prefactor_r[j, :]

    # Peierls inter-orbital displacement (tau_L - tau_Lp in Cartesian)
    # Matches lwann exactly: distances = plist[:,None,:] - plist[None,:,:]
    if peierlscorrection:
        distances   = self.plist[:, None, :] - self.plist[None, :, :]
        ri_minus_rj = np.einsum('id,abi->abd', self.rvec, distances)  # (norb,norb,3)

    # ── spin loop ─────────────────────────────────────────────────────────
    for ispin in range(self.spins):
        hr = self.hrlist[ispin]     # (nrpts, norb, norb) complex

        # ── temperature loop ──────────────────────────────────────────────
        for iT, sp in enumerate(se_params_list):
            logger.info("  T={:.2f} K  spin {}/{}  iT {}/{}".format(
                sp.temperature, ispin+1, self.spins, iT+1, nt))

            # unpack k-independent QP parameters for this (T, spin)
            Zinv_orb     = sp.Zinv_orb                          # (norb, norb) real
            Gamma_orb    = sp.Gamma_orb                          # (norb, norb) real, symmetric (may be off-diagonal-aware)
            mu_dmft      = sp.mu

            # Sigma0dc_orb = Sigma0_orb (full matrix, off-diagonal-aware)
            #              + dc_orb on the diagonal only (w2dynamics double
            #              counting is always diagonal in orbital space:
            #              shape (norb,) per ineq. site, never (norb,norb)).
            Sigma0dc_orb = sp.Sigma0_orb.copy()                  # (norb, norb) complex
            diag_idx_T   = np.arange(self.nproj)
            Sigma0dc_orb[diag_idx_T, diag_idx_T] += sp.dc_orb

            # off-diagonal scattering diagnostic accumulator (Approximation 3)
            offdiag_diags = []   # list of OffdiagDiag, one per k (see _per_k_dmft)

            # ── k loop ────────────────────────────────────────────────────
            for ik in range(self.nkp):
                kpt   = self.kpoints[ik]                       # (3,) fractional
                rdotk = 2.0 * np.pi * (self.rlist @ kpt)      # (nrpts,)
                ee    = np.exp(1j * rdotk) / self.rmultiplicity # (nrpts,)

                # FT: H(k)  shape (norb, norb)
                hk_orb = np.einsum('r,rij->ij', ee, hr)

                # FT: velocity H(k)  shape (norb, norb, 3)
                hvk_orb = np.einsum('dr,r,rij->ijd', 1j * prefactor_r, ee, hr)

                # Peierls correction (matching lwann): -1j * (tau_L - tau_Lp) * H(k)
                if peierlscorrection:
                    hvk_orb += -1j * hk_orb[:, :, None] * ri_minus_rj

                # FT: curvature H(k)  shape (norb, norb, 6)
                hck_orb = np.einsum('dr,r,rij->ijd', -prefactor_r2, ee, hr)

                # per-k DMFT renormalisation
                (eps_bare, E_QP, gamma_QP, Z_QP, bshift,
                 vel2diag, Bvel2diag, vel2, Bvel2, diag) = _per_k_dmft(
                    hk_orb, hvk_orb, hck_orb,
                    Zinv_orb, Sigma0dc_orb, Gamma_orb,
                    mu_dmft, sp.beta, self.ortho)

                offdiag_diags.append(diag)

                # store spectral QP arrays.
                #
                # IMPORTANT: gamma_QP returned by _per_k_dmft is already the
                # fully Z-renormalised physical rate (tex sec:QPdiag); but
                # LinReTraCe's own Fortran core (src_linretrace/input.f90)
                # performs `sct%gam = sct%gam * sct%zqp` on EVERY scattering
                # HDF5 it reads, regardless of who wrote it or how. Writing
                # gamma_QP unmodified to 'scatrate' would let LinReTraCe
                # multiply it by Z_QP a SECOND time downstream, silently
                # double-suppressing the scattering rate (Z^2*Gamma instead
                # of Z*Gamma in the single-orbital case). Dividing by Z_QP
                # here exactly compensates that read-in convention, so the
                # rate LinReTraCe actually uses downstream is gamma_QP, as
                # intended. Z_QP itself is always > 0 (it is a sum of squared
                # eigenvector components, Eq. Ztilde_final), so this division
                # is safe.
                gamma_out    [iT, ispin, ik, :] = gamma_QP / Z_QP
                Z_out        [iT, ispin, ik, :] = Z_QP
                bandshift_out[iT, ispin, ik, :] = bshift

            # off-diagonal scattering diagnostic for this (T, spin),
            # see tex sec:validity item 5 and structure/checkfit.py
            _report_offdiag_scattering_diagnostic(
                offdiag_diags, self.kpoints, sp.temperature, ispin, self.spins)

    # ── floor on scattering rates (FullScattering warns below 1e-14) ──────
    # Floor applies to the value about to be WRITTEN to 'scatrate'
    # (gamma_QP / Z_QP), consistent with FullScattering's own internal
    # check on self.scattering, and with what LinReTraCe reads in before
    # its own *zqp multiplication.
    n_floor = int(np.sum(gamma_out < 1e-14))
    if n_floor > 0:
        logger.warning(
            "{} scattering rate entries below 1e-14 eV; "
            "applying floor.".format(n_floor))
        gamma_out = np.maximum(gamma_out, 1e-14)

    # ── write scattering HDF5 ─────────────────────────────────────────────
    scatobj.defineScatteringRates(
        gamma_out.astype(np.float64),
        qpweight  = Z_out.astype(np.float64),
        bandshift = bandshift_out.astype(np.float64),
    )
    # FullScattering.createOutput() prompts interactively ("Overwrite?
    # (y/N)") if `output` already exists, and just returns (None, exactly
    # like on success -- there is no way to tell the two apart from its
    # return value) without writing anything if the user declines.  It
    # does already log either outcome correctly itself ("Successfully
    # created scattering file: ..." or "File will not be overwritten.");
    # what was missing is that the line below used to fire unconditionally
    # afterwards regardless of which happened.  Detect an actual write,
    # without touching scattering/fullscattering.py (kept unmodified
    # throughout this interface, see tex sec:impl), by comparing the
    # file's mtime before and after the call -- a declined prompt leaves
    # the pre-existing file's mtime untouched, a real write always
    # advances it (including the "file did not exist before" case, where
    # mtime_before is None).
    mtime_before = os.path.getmtime(output) if os.path.isfile(output) else None
    scatobj.createOutput(output)
    written = os.path.isfile(output) and (
        mtime_before is None or os.path.getmtime(output) > mtime_before)
    if written:
        logger.info("computeDMFT: written -> {}".format(output))
    else:
        logger.warning(
            "computeDMFT: scattering file NOT written -- the overwrite "
            "prompt above was likely declined: {}".format(output))
