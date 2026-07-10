#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
structure/wannier_disorder.py
=============================
Ab-initio DISORDER input for LinReTraCe from a two-point disorder archive
(deltaH.h5 of the disorder-Wannier pipeline).  Grafted onto
`Wannier90Calculation` in lwann, exactly like the DMFT extension:

    from structure.wannier_disorder import computeDisorder as _computeDisorder
    Wannier90Calculation.computeDisorder = _computeDisorder

The LinReTraCe Fortran core is NOT touched.  Outputs:

  1. a scattering HDF5 (FullScattering layout) with momentum- and
     band-resolved on-shell disorder rates, virtual-crystal (Sigma^(1))
     bandshifts, and Z = 1 (elastic disorder is frequency-independent);
  2. optionally a DRESSED energy HDF5: a copy of the plain-lwann energy
     file whose transition matrix elements carry the static disorder
     ladder (diffuson) vertex, M~^{d1 d2}_{nm} = Re[conj(v^{d1}_{nm})
     Lambda^{d2}_{nm}] -- only ONE side is dressed (both sides would
     double count the ladder);
  3. optionally a LAMBDA HDF5 with the vertex enhancements themselves
     (diagonal ratios, complex vertex, inter-band pairs) for analysing
     their momentum and band structure.

Physics conventions
-------------------
* Rates (engine='born'|'scba'|'ata'): band-basis on-shell
      Gamma_kn = sum_{p,c,m} w_c |dV^band_c(kn -> k'm)|^2
                 Z_m Gt_{k'm} / [(a_kn - a_{k'm})^2 + Gt_{k'm}^2],
  k' = k - q_p, with Gt the total renormalised rail width (iterated to
  self-consistency for 'scba'; 'ata' resums per-configuration multiple
  scattering on eta rails and captures impurity bound-state trends).
  Elastic => T-independent; replicated over the requested temperature
  axis, unless combined with --dmft (see below).
* Vertex evaluation energy (vertex_eval flag): the weight multiplying
  lambda_k(eps) in the transport kernel is (-f')(eps) |G_k(eps)|^2, the
  product of the thermal window (width T at mu) and the state's own
  Lorentzian (width Gt at a_kn).  Because (-f') decays exponentially
  while the Lorentzian tail decays only polynomially, this product is in
  general BIMODAL for |a_kn| >> T: a tail lobe at eps ~ 0 (height
  ~1/a^2) and an on-shell lobe at eps ~ a_kn (height ~exp(-a/T)/Gt^2);
  the on-shell lobe dominates iff a/T < 2 ln(a/Gt).
    - 'fermi'  : lambda solved once per orbit at eps = 0.  Correct for
      metals (contributing states have |a| <~ max(T, Gt), both choices
      then agree to leading order) and for the LOW-T tail-dominated
      (resistivity-saturation) regime of narrow-gap semiconductors.
    - 'onshell': lambda solved at eps = a_kn per state (pair midpoint
      for inter-band elements), energy-binned for cost.  Correct in the
      ACTIVATED regime of semiconductors, where carriers live at the
      band edges and eps = 0 lies in the gap where no state exists;
      also retains the particle-hole asymmetry of the vertex that a
      single lambda(0) suppresses (relevant for L12/Seebeck).
  The DISAGREEMENT between the two settings is itself the diagnostic
  for the bimodal crossover regime, where any single static dressing is
  strained; there, use the full disorder_transport machinery instead.
* Conserved-mode deflation: for elastic scattering the per-energy Ward
  identity places a near-unit eigenvalue of the RA ladder kernel (the
  conserved-density/diffusion mode) at EVERY energy.  The current
  source is orthogonal to it by parity for inversion-symmetric,
  self-averaging correlators.  The orbit solver therefore uses a
  truncated-SVD pseudo-inverse with a RELATIVE cut (default 1e-2): the
  near-null direction is dropped, and the source amplitude that would
  have entered the dropped directions ("leak") is monitored -- a large
  leak means broken parity protection (non-self-averaging ensemble) and
  is reported as a warning.
* Fortran read-in convention (as in the DMFT interface): the core
  multiplies scatrate by qpweight on read.  Disorder-only: Z = 1 and
  rates are written as computed.  Combined with --dmft: written rate =
  gamma_dmft_QP / Z_QP + Gamma_dis, so the rate used downstream is
  Z*(Gamma0_dmft + Gamma_dis).

Combination with --dmft
-----------------------
Static Hermitian parts add BEFORE the rotation:
H~(k) = H(k) + Sigma0_dmft + dc + Sigma^(1)_dis(k) - mu_dmft, solved as
the same generalised eigenvalue problem with Zinv_dmft as in
wannier_dmft (Sigma^(1) is k-dependent, hence added per k).  The rails
of the disorder rate and of the ladder carry the TOTAL renormalised
width and the spectral weight Z (coherent part only, consistent with
LinReTraCe); the rung remains disorder-only (no e-e vertex by
assumption), so lambda acquires a parametric T-dependence through the
DMFT rates and is computed per temperature.  The DMFT temperature axis
governs.
"""

import os
import sys
import shutil
import logging
from collections import defaultdict

import numpy as np
import h5py

logger = logging.getLogger(__name__)

from scattering.fullscattering import FullScattering

_KB_EV = 8.61733034e-5


# ---------------------------------------------------------------------------
def _import_dwannier(path_hint=None):
    """Import the disorder pipeline modules (disorder_selfenergy, dwannier,
    disorder_transport)."""
    for p in (path_hint, os.environ.get("DWANNIER_PATH")):
        if p and p not in sys.path:
            sys.path.insert(0, p)
    try:
        import disorder_selfenergy as ds
        import disorder_transport as dtmod
        return ds, dtmod
    except ImportError as e:
        raise ImportError(
            "computeDisorder needs the disorder-Wannier pipeline modules "
            "(disorder_selfenergy.py, disorder_transport.py, dwannier.py). "
            "Provide their location via --dwannier-path or DWANNIER_PATH. "
            "({})".format(e))


def _mesh_index_map(kpoints, shape):
    """LRT k-point index -> uniform-grid index; requires an unshifted
    Monkhorst-Pack grid."""
    n = np.asarray(shape, int)
    x = np.mod(np.asarray(kpoints, float), 1.0) * n[None, :]
    if np.abs(x - np.rint(x)).max() > 1e-6:
        raise ValueError(
            "k-mesh is not a uniform unshifted {}x{}x{} grid; --disorder "
            "requires a Gamma-centred Monkhorst-Pack mesh.".format(*n))
    ix = np.rint(x).astype(int) % n
    return ix[:, 0] * n[1] * n[2] + ix[:, 1] * n[2] + ix[:, 2]


def _orbits_from_jshift(nkp, jshift):
    """Partition k indices into closed q_p-transfer orbits."""
    seen = np.full(nkp, -1, int)
    orbits = []
    for ik in range(nkp):
        if seen[ik] >= 0:
            continue
        orb, stack = [], [ik]
        while stack:
            i = stack.pop()
            if seen[i] >= 0:
                continue
            seen[i] = len(orbits)
            orb.append(i)
            stack.extend(int(j) for j in jshift[i])
        orbits.append(sorted(orb))
    return orbits


# ---------------------------------------------------------------------------
def _onshell_rates(a, jshift, V2, w_norm, Zw, engine='scba',
                   gamma_other=None, gamma_seed=1e-2, max_iter=400,
                   tol=1e-12, rate_eval='onshell'):
    """On-shell Born/SCBA disorder rates (bare, i.e. before the Fortran
    core's own *Z).  a, Zw, gamma_other: (nkp, nb).  V2: (nkp, N, nb, nb)
    ensemble-averaged |dV^band|^2 with V2[ik, p, n, m] coupling (ik, n) to
    (jshift[ik,p], m)."""
    nkp, nb = a.shape
    g_oth = np.zeros_like(a) if gamma_other is None else gamma_other
    G = np.full_like(a, gamma_seed)
    mix = 0.5
    res = np.inf
    for it in range(max_iter if engine == 'scba' else 1):
        rail = Zw * (G + g_oth)
        Gn = np.zeros_like(G)
        for p in range(jshift.shape[1]):
            j = jshift[:, p]
            num = (Zw * rail)[j]                       # (nkp, nb) at k'
            aev = a if rate_eval == 'onshell' else np.zeros_like(a)
            den = (aev[:, :, None] - a[j][:, None, :]) ** 2 \
                + rail[j][:, None, :] ** 2             # (nkp, nb, nb)
            Gn += np.einsum('knm,knm->kn', V2[:, p], num[:, None, :] / den)
        if engine != 'scba':
            return np.maximum(Gn, 0.0)
        res_new = float(np.abs(Gn - G).max())
        if res_new < tol:
            return Gn
        if res_new > res:
            mix = max(0.1, 0.6 * mix)
        G = (1 - mix) * G + mix * Gn
        res = res_new
    logger.warning("    SCBA on-shell rates: residual %.1e after %d "
                   "iterations (small/non-self-averaging ensembles may "
                   "not admit a stable fixed point)", res, max_iter)
    return G


def _ata_rates(a, orbits, jshift, dVband, w_conf, Zw, gamma_other=None,
               eta=None):
    """Average-T-matrix on-shell rates: per configuration, multiple
    scattering resummed within each q_p orbit on Lorentzian rails.
    Captures bound-state trends beyond second order."""
    nkp, nb = a.shape
    g_oth = np.zeros_like(a) if gamma_other is None else gamma_other
    born = _onshell_rates(a, jshift, _v2_from_dvband(dVband, w_conf,
                                                     jshift.shape[1], nb),
                          w_conf, Zw, engine='born', gamma_other=g_oth)
    rail = np.maximum(Zw * (born + g_oth), 1e-4) if eta is None else \
        np.full_like(a, eta)
    G = np.zeros((nkp, nb))
    for orb in orbits:
        no = len(orb)
        pos = {ik: i for i, ik in enumerate(orb)}
        for ik in orb:
            for n in range(nb):
                e = a[ik, n]
                gd = np.concatenate(
                    [Zw[jk] / (e - a[jk] + 1j * rail[jk]) for jk in orb])
                for c, wc in enumerate(w_conf):
                    Vb = np.zeros((no * nb, no * nb), complex)
                    for i, jk in enumerate(orb):
                        for p in range(jshift.shape[1]):
                            b = pos[int(jshift[jk, p])]
                            Vb[i*nb:(i+1)*nb, b*nb:(b+1)*nb] += \
                                dVband[jk][p][c]
                    T = np.linalg.solve(
                        np.eye(no * nb, dtype=complex) - Vb * gd[None, :],
                        Vb)
                    i0 = pos[ik]
                    G[ik, n] += -wc * np.imag(T[i0*nb + n, i0*nb + n])
    return np.maximum(G, 0.0)


def _v2_from_dvband(dVband, w_conf, N, nb):
    nkp = len(dVband)
    V2 = np.zeros((nkp, N, nb, nb))
    for ik in range(nkp):
        for p in range(N):
            for c, wc in enumerate(w_conf):
                V2[ik, p] += wc * np.abs(dVband[ik][p][c]) ** 2
    return V2


# ---------------------------------------------------------------------------
def _solve_vertex_orbit(orb, jshift, dVband, w_conf, a, Gt, Zw, vband,
                        ebar, svd_cut, stats):
    """Static disorder ladder on ONE q_p orbit at energy ebar, full band-
    matrix form (RA section):

        Lambda_{nm}(k) = v_{nm}(k) + sum_{p,c} w_c [ dV_c(k,k')
                          G^R(k') Lambda(k') G^A(k') dV_c(k,k')^dag ]_{nm}

    with band-diagonal rails G^{R/A}_m(k') = Z_m/(ebar - a_{k'm} +/-
    i Gt_{k'm}).  Solved for the three Cartesian sources via a
    truncated-SVD pseudo-inverse of (1-K) with RELATIVE cut svd_cut:
    the near-unit conserved-density (diffusion) direction implied by the
    per-energy Ward identity is dropped; `stats` accumulates the number
    of dropped directions, the worst source leak into them, and the
    effective condition number.  Returns (3, len(orb), nb, nb) complex.
    """
    no = len(orb)
    nb = a.shape[1]
    dim = no * nb * nb
    pos = {ik: i for i, ik in enumerate(orb)}
    GR = {jk: Zw[jk] / (ebar - a[jk] + 1j * Gt[jk]) for jk in orb}
    K = np.zeros((dim, dim), complex)
    for i, ik in enumerate(orb):
        for p in range(jshift.shape[1]):
            jk = int(jshift[ik, p])
            b = pos[jk]
            gr = GR[jk]
            blk = np.zeros((nb * nb, nb * nb), complex)
            for c, wc in enumerate(w_conf):
                V = dVband[ik][p][c]
                A = V * gr[None, :]                    # V @ diag(G^R)
                # B = diag(G^A) V^dag ; row-major vec(AXB) = (A kron B^T) x
                # with B^T[m, m'] = conj(V[m, m']) conj(gr[m'])
                Bt = np.conj(V) * np.conj(gr)[None, :]
                blk += wc * np.kron(A, Bt)             # row-major vec(A L B)
            K[i*nb*nb:(i+1)*nb*nb, b*nb*nb:(b+1)*nb*nb] += blk
    M = np.eye(dim, dtype=complex) - K
    U, S, Vh = np.linalg.svd(M)
    keep = S > svd_cut * S.max()
    stats["dropped"] += int(np.sum(~keep))
    stats["cond"] = max(stats["cond"], float(S.max() / S[keep].min()))
    out = np.empty((3, no, nb, nb), complex)
    for al in range(3):
        src = np.stack([vband[ik][:, :, al] for ik in orb]).reshape(dim)
        coef = U.conj().T @ src
        if (~keep).any():
            leak = float(np.abs(coef[~keep] / S[~keep]).max())
        x = Vh.conj().T[:, keep] @ (coef[keep] / S[keep])
        if (~keep).any():
            stats["leak"] = max(stats["leak"],
                                leak / max(np.linalg.norm(x), 1e-300))
        out[al] = x.reshape(no, nb, nb)
    return out



def _cooperon_orbitQ(orb, Qv, kpoints, shape, jshift, dVband, w_conf,
                     a, Gt, Zw, vband, ebar, mesh_index):
    """Maximally crossed (n >= 2) additive current vertex on one q_p orbit
    at fixed total momentum Q and energy ebar (static band-diagonal RA
    rails).  Transplant of the machine-validated crossed-fan transfer
    operator of disorder_transport._crossed_fan_Q:

        O[(u),(u-q_p)] = sum_c w_c [dV_c(u,u') G^R(u')] kron
                                   [dV_c(Q-u,(Q-u)+q_p) G^A((Q-u)+q_p)],
        dV_c(Q-u,(Q-u)+q_p) = dV_c((Q-u)+q_p -> Q-u)^dagger,

    fan sum n>=2 = O^2 (1-O)^{-1}, closed per k in the (k, Q-k) block with
    the twisted boundary.  Returns {ik: C (3, nb, nb)} with the additive
    vertex C^beta(k) defined by  t_k = Tr[v^alpha G^R C^beta G^A], i.e.
    C_{il} = P_{il}/G^A_l(k),
    P_{il} = sum_{j k'} G^A_j(Q-k) X_{(i j),(k' l)} v^beta(Q-k)_{k' j}.
    Momenta outside the orbit contribute nothing (no closure)."""
    no = len(orb)
    nb = a.shape[1]
    d = nb * nb
    pos = {ik: i for i, ik in enumerate(orb)}
    GR = {jk: Zw[jk] / (ebar - a[jk] + 1j * Gt[jk]) for jk in orb}
    iQmk = {}
    closes = False
    for ik in orb:
        j = mesh_index(Qv - kpoints[ik])
        iQmk[ik] = j
        if j in pos:
            closes = True
    if not closes:
        return {}
    # G^A at arbitrary mesh points (lower rail lives at Q-u, may leave orb)
    def GA_at(jk):
        return np.conj(Zw[jk] / (ebar - a[jk] + 1j * Gt[jk]))
    O = np.zeros((no * d, no * d), complex)
    for iu in orb:
        i = pos[iu]
        for p in range(jshift.shape[1]):
            iup = int(jshift[iu, p])
            if iup not in pos:
                continue
            b = pos[iup]
            il = iQmk[iu]
            ilp = mesh_index(Qv - kpoints[iup])       # (Q-u) + q_p
            gr = GR[iup]
            ga = GA_at(ilp)
            blk = np.zeros((d, d), complex)
            for c, wc in enumerate(w_conf):
                A = dVband[iu][p][c] * gr[None, :]
                B = dVband[ilp][p][c].conj().T * ga[None, :]
                blk += wc * np.kron(A, B)
            O[i*d:(i+1)*d, b*d:(b+1)*d] += blk
    X = O @ O @ np.linalg.solve(np.eye(no * d, dtype=complex) - O,
                                np.eye(no * d, dtype=complex))
    out = {}
    for ik in orb:
        jf = iQmk[ik]
        if jf not in pos:
            continue
        Xb = X[pos[ik]*d:(pos[ik]+1)*d,
               pos[jf]*d:(pos[jf]+1)*d].reshape(nb, nb, nb, nb)
        gaQ = GA_at(jf)                               # (nb,) at Q-k
        gak = np.conj(GR[ik])                         # (nb,) at k
        C = np.empty((3, nb, nb), complex)
        for be in range(3):
            vbQ = vband[jf][:, :, be]
            # P_il = sum_{j,kp} gaQ_j * Xb[i,j,kp,l] * vbQ[kp,j]
            P = np.einsum('j,ijkl,kj->il', gaQ, Xb, vbQ)
            C[be] = P / gak[None, :]
        out[ik] = C
    return out


# ---------------------------------------------------------------------------
def computeDisorder(self,
                    archive_path,
                    energy_hdf5,
                    scat_output='lrtc-disorder-scat.hdf5',
                    temps=(100.0, 500.0, 5, False),
                    engine='scba',
                    vca=True,
                    vertex=True,
                    vertex_eval='onshell',
                    rate_eval='onshell',
                    vertex_channels='ladder',
                    gamma_add=(0.0, 0.0),
                    dressed_per_t=False,
                    dressed_output='lrtc-wann-dressed.hdf5',
                    lambda_output=None,
                    coherence_window=6.0,
                    svd_cut=1e-2,
                    dwannier_path=None,
                    se_params_list=None,
                    peierlscorrection=True,
                    interactive=True):
    """Generate LinReTraCe disorder input from a deltaH.h5 archive.  See
    the module docstring for physics and lwann --help for the CLI mapping.

    vertex_eval : 'onshell' (lambda at eps = a_kn / pair midpoints,
                  energy-binned) or 'fermi' (lambda at eps = 0, one solve
                  per orbit).  See module docstring for the regime map.
    rate_eval   : 'onshell' (Gamma_kn = -Im Sigma at eps = a_kn, the
                  quasiparticle width) or 'fermi' (-Im Sigma at eps = 0,
                  governing the tail spectral weight at mu).  Keep matched
                  with vertex_eval for a consistent linearised pair.
    vertex_channels : 'ladder' (diffuson only, default), 'cooperon'
                  (maximally crossed only), 'both'.  The Cooperon requires
                  an inelastic dephasing rate on the rails (gamma_add or
                  DMFT); see the notes for the cutoff physics.
    gamma_add   : (c0, c2) phenomenological INELASTIC rate
                  Gamma_inel(T) = c0 + c2*T^2  [eV, eV/K^2] added to the
                  rails of rates and vertex AND to the written scatrate
                  (per temperature step).  This keeps rates, ladder,
                  Cooperon cutoff and kernels mutually consistent --
                  preferred over a downstream ScatteringOffset when
                  vertex corrections are on.  Makes the vertex (and hence
                  the dressing) temperature dependent.
    dressed_per_t : when the dressing is T-dependent (combined mode or
                  gamma_add != 0), write one dressed energy file per
                  temperature ('<name>_T####.hdf5') instead of only the
                  first temperature.
    svd_cut     : RELATIVE truncated-SVD cut deflating the near-unit
                  conserved-density ladder direction (default 1e-2).
    """
    ds, dtmod = _import_dwannier(dwannier_path)
    if vertex_eval not in ('onshell', 'fermi'):
        raise ValueError("vertex_eval must be 'onshell' or 'fermi'")
    if rate_eval not in ('onshell', 'fermi'):
        raise ValueError("rate_eval must be 'onshell' or 'fermi'")
    if vertex_channels not in ('ladder', 'cooperon', 'both'):
        raise ValueError("vertex_channels: 'ladder'|'cooperon'|'both'")
    gamma_add = (float(gamma_add[0]), float(gamma_add[1]))
    if rate_eval != vertex_eval and vertex:
        logger.warning("rate_eval=%s and vertex_eval=%s differ: rails and"
                       " ladder are then not a matched linearised pair;"
                       " keep them equal unless you have a reason.",
                       rate_eval, vertex_eval)

    logger.info("=" * 64)
    logger.info("computeDisorder: ab-initio disorder -> LinReTraCe")
    logger.info("  archive     : {}".format(archive_path))
    logger.info("  energy file : {}".format(energy_hdf5))
    logger.info("  engine={} vca={} vertex={} vertex_eval={}".format(
        engine, vca, vertex, vertex_eval))

    arch = ds.DisorderArchive(archive_path)
    wc = np.asarray(arch.w, float)
    nc = len(wc)
    ref = arch.reference(average_potential=False)
    norb = ref.H0(np.zeros(3)).shape[0]
    if norb != self.nproj:
        raise ValueError("archive n_p={} != Wannier90 nproj={}".format(
            norb, self.nproj))

    # ── mesh commensurability and index maps ──────────────────────────────
    shape = (self.nkx, self.nky, self.nkz)
    Qf = np.asarray(arch.Q_frac, float)
    qint = Qf * np.asarray(shape, float)[None, :]
    if np.abs(qint - np.rint(qint)).max() > 1e-8:
        raise ValueError(
            "k-mesh {} is not commensurate with the archive transfers "
            "q_p = {}; each q_p * (nkx,nky,nkz) must be integer "
            "(choose --kmesh accordingly).".format(
                shape, np.round(Qf, 4).tolist()))
    if getattr(self, "kshift", False):
        raise ValueError("--disorder requires an unshifted k-mesh.")
    lin = _mesh_index_map(self.kpoints, shape)         # lrt ik -> grid id
    setup = dtmod.TransportSetup(arch, ref, shape)     # grid-ordered caches
    grid2lrt = np.empty(self.nkp, int)
    grid2lrt[lin] = np.arange(self.nkp)
    jshift = grid2lrt[setup.jshift[lin]]               # lrt-indexed
    N = setup.N
    orbits = _orbits_from_jshift(self.nkp, jshift)
    ineg = grid2lrt[np.array([ds._mesh_index(-self.kpoints[i], shape)
                              for i in range(self.nkp)])]
    logger.info("  {} configurations, {} transfers, {} orbits "
                "(max size {})".format(nc, N, len(orbits),
                                       max(len(o) for o in orbits)))

    # ── energy file / FullScattering / temperature axis ───────────────────
    scatobj = FullScattering(energy_hdf5)
    if not interactive:
        scatobj.inputmethod = lambda *_a, **_k: 'y'
    spins, nkp_e, nbands = scatobj.getDependencies()
    if nkp_e != self.nkp or nbands != norb or spins != self.spins:
        raise ValueError(
            "energy file (nkp={}, nbands={}, spins={}) inconsistent with "
            "this lwann run (nkp={}, nproj={}, spins={}); regenerate it "
            "on the same --kmesh.".format(nkp_e, nbands, spins, self.nkp,
                                          norb, self.spins))
    mu = float(scatobj.mudft)
    combined = bool(se_params_list)
    if combined:
        Ts = np.array([sp.temperature for sp in se_params_list], float)
    else:
        tmin, tmax, ntt, tlog = temps
        Ts = (np.logspace(np.log10(tmin), np.log10(tmax), int(ntt))
              if tlog else np.linspace(tmin, tmax, int(ntt)))
    nt = len(Ts)
    scatobj.defineTemperatures(float(Ts[0]), float(Ts[-1]), nt, tlog=False)
    scatobj.temps = Ts
    scatobj.betas = 1.0 / (Ts * _KB_EV)
    scatobj.mus = np.full(nt, mu)
    scatobj.tmin, scatobj.tmax = float(Ts[0]), float(Ts[-1])

    # ── FT prefactors (identical to computeDMFT / computeHamiltonian) ─────
    prefactor_r = np.einsum('id,ri->dr', self.rvec, self.rlist)
    ri_minus_rj = None
    if peierlscorrection:
        distances = self.plist[:, None, :] - self.plist[None, :, :]
        ri_minus_rj = np.einsum('id,abi->abd', self.rvec, distances)

    dir1 = [0, 1, 2, 0, 0, 1]
    dir2 = [0, 1, 2, 1, 2, 2]
    ndir = 3 if self.ortho else 9

    gamma_out = np.zeros((nt, spins, self.nkp, norb))
    Z_out = np.ones((nt, spins, self.nkp, norb))
    bshift_out = np.zeros((nt, spins, self.nkp, norb))
    lamD = np.ones((max(nt, 1), spins, self.nkp, norb, ndir))
    vertD = np.zeros((max(nt, 1), spins, self.nkp, norb, 3), complex)
    MdiagD = np.zeros((max(nt, 1), spins, self.nkp, norb, ndir))
    inter = defaultdict(dict)       # (iT,ispin,ik) -> {(n,m): (Lam3, M6)}
    guard_skipped = 0

    if combined:
        import scipy.linalg as sla

    per_t = combined or (gamma_add[0] != 0.0 or gamma_add[1] != 0.0)
    ladder_on = vertex and vertex_channels in ('ladder', 'both')
    coop_on = vertex and vertex_channels in ('cooperon', 'both')
    tloop = range(nt) if per_t else [0]
    for ispin in range(self.spins):
        hr = self.hrlist[ispin]
        for iT in tloop:
            sp = se_params_list[iT] if combined else None
            # ── rotation, band energies, band-basis velocities ────────────
            a = np.zeros((self.nkp, norb))              # E - mu
            eps_bare = np.zeros((self.nkp, norb))
            Zw = np.ones((self.nkp, norb))
            g_dmft = np.zeros((self.nkp, norb))
            Ulist = [None] * self.nkp
            vband = [None] * self.nkp
            g_phen = gamma_add[0] + gamma_add[1] * float(Ts[iT]) ** 2
            if combined:
                S0dc = sp.Sigma0_orb.copy()
                S0dc[np.arange(norb), np.arange(norb)] += sp.dc_orb
            h0dev = 0.0
            for ik in range(self.nkp):
                kpt = self.kpoints[ik]
                rdotk = 2.0 * np.pi * (self.rlist @ kpt)
                ee = np.exp(1j * rdotk) / self.rmultiplicity
                hk = np.einsum('r,rij->ij', ee, hr)
                hv = np.einsum('dr,r,rij->ijd', 1j * prefactor_r, ee, hr)
                if peierlscorrection:
                    hv += -1j * hk[:, :, None] * ri_minus_rj
                if ik % max(1, self.nkp // 5) == 0:
                    h0dev = max(h0dev,
                                float(np.abs(hk - ref.H0(kpt)).max()))
                eps_bare[ik] = np.linalg.eigvalsh(hk)
                S1 = arch.sigma1(kpt) if vca else np.zeros((norb, norb))
                if combined:
                    Ht = hk + S0dc + S1 - sp.mu * np.eye(norb)
                    Ht = 0.5 * (Ht + Ht.conj().T)
                    E, W = sla.eigh(Ht, sp.Zinv_orb)
                    Zw[ik] = np.real(np.einsum('an,an->n', W.conj(), W))
                    g_dmft[ik] = np.real(np.einsum(
                        'an,ab,bn->n', W.conj(), sp.Gamma_orb, W)) / Zw[ik]
                    a[ik] = E                            # mu-referenced
                    bshift_out[iT, ispin, ik] = (E / Zw[ik] - eps_bare[ik]
                                                 + sp.mu)
                else:
                    Hs = hk + S1
                    E, W = np.linalg.eigh(0.5 * (Hs + Hs.conj().T))
                    a[ik] = E - mu
                    bshift_out[:, ispin, ik] = E - eps_bare[ik]
                Ulist[ik] = W
                vband[ik] = np.einsum('an,abd,bm->nmd', W.conj(), hv, W)
            if ispin == 0 and iT == tloop[0]:
                msg = "  archive/lwann H0(k) consistency: {:.2e} eV".format(
                    h0dev)
                (logger.warning if h0dev > 1e-6 else logger.info)(msg)

            # ── band-basis fluctuations and correlator diagnostics ────────
            dVband = [[[Ulist[ik].conj().T @ setup.dV[lin[ik], p, c]
                        @ Ulist[int(jshift[ik, p])]
                        for c in range(nc)] for p in range(N)]
                      for ik in range(self.nkp)]
            V2 = _v2_from_dvband(dVband, wc, N, norb)
            V2_export, a_export, Gt_export = V2, None, None
            Wdiag = np.einsum('kpnn->kn', V2)
            asym = float(np.abs(Wdiag - Wdiag[ineg]).max()
                         / max(Wdiag.max(), 1e-300))
            (logger.warning if asym > 0.05 else logger.info)(
                "  correlator parity asymmetry {:.1%}{}".format(
                    asym, " > 5%: ensemble not self-averaging; parity "
                    "protection of the current channel is weakened -- "
                    "treat lambda qualitatively and monitor the leak "
                    "diagnostic" if asym > 0.05 else
                    " (parity protection OK)"))

            # ── on-shell rates ────────────────────────────────────────────
            g_oth = g_dmft + g_phen
            if engine == 'ata':
                Gdis = _ata_rates(a, orbits, jshift, dVband, wc, Zw,
                                  gamma_other=g_oth)
            else:
                Gdis = _onshell_rates(a, jshift, V2, wc, Zw, engine,
                                      gamma_other=g_oth,
                                      rate_eval=rate_eval)
            Gt = Zw * (Gdis + g_oth)
            a_export, Gt_export = a.copy(), Gt.copy()
            lvl = float(np.median(np.abs(np.diff(np.sort(a.ravel())))))
            if float(np.median(Gt)) < lvl:
                logger.warning(
                    "  median rail width {:.3g} eV < median level spacing "
                    "{:.3g} eV: on-shell quantities are k-mesh limited; "
                    "densify --kmesh.".format(float(np.median(Gt)), lvl))
            logger.info("  Gamma_dis in [{:.4g}, {:.4g}] eV "
                        "(median {:.4g})".format(Gdis.min(), Gdis.max(),
                                                 float(np.median(Gdis))))
            tsl = [iT] if per_t else list(range(nt))
            for jt in tsl:
                gamma_out[jt, ispin] = Gdis + g_oth
                Z_out[jt, ispin] = Zw

            # ── static vertex: diffuson ladder and/or Cooperon ───────────
            if not vertex:
                continue
            stats = {"dropped": 0, "leak": 0.0, "cond": 1.0, "pairs": 0}
            iTl = iT if per_t else 0
            ebin = max(float(Gt.min()) / 3.0, 1e-4)
            # target (k, n, m) elements grouped per orbit and energy bin
            orb_targets = []
            for orb in orbits:
                targ = defaultdict(list)
                for ik in orb:
                    for n in range(norb):
                        e_n = 0.0 if vertex_eval == 'fermi' else a[ik, n]
                        targ[int(round(e_n / ebin))].append((ik, n, n))
                        for m in range(n + 1, norb):
                            if (abs(a[ik, n] - a[ik, m]) < coherence_window
                                    * (Gt[ik, n] + Gt[ik, m])):
                                e_p = 0.0 if vertex_eval == 'fermi' else \
                                    0.5 * (a[ik, n] + a[ik, m])
                                targ[int(round(e_p / ebin))].append(
                                    (ik, n, m))
                                stats["pairs"] += 1
                orb_targets.append((orb, dict(targ)))

            # Cooperon channel: physical prechecks + Q-parallel pass
            Cmap = {}
            if coop_on:
                g_inel_med = float(np.median(Zw * g_oth))
                if g_inel_med < 1e-6:
                    raise ValueError(
                        "vertex_channels includes the Cooperon but no "
                        "inelastic dephasing rate is on the rails (median "
                        "{:.1e} eV): the Q -> 0 pole is unregularised. "
                        "Supply --gamma-add C0 C2 (or --dmft); elastic "
                        "disorder does not dephase the Cooperon."
                        .format(g_inel_med))
                # length scales: mean free path, dephasing length,
                # supercell size, Q-mesh resolution
                Bmat = 2.0 * np.pi * np.linalg.inv(self.rvec).T
                vd0 = np.stack(
                    [np.real(vband[ik][np.arange(norb), np.arange(norb)])
                     for ik in range(self.nkp)])         # (nk, nb, 3)
                win = np.abs(a) < max(float(np.median(Gt)), 0.05)
                v2all = np.sum(vd0 ** 2, axis=-1)
                vF2 = float(np.mean(v2all[win])) if win.any() \
                    else float(np.mean(v2all))
                g_el_med = max(float(np.median(Zw * Gdis)), 1e-6)
                tau_el = 1.0 / (2.0 * g_el_med)
                tau_phi = 1.0 / (2.0 * g_inel_med)
                ell = np.sqrt(vF2) * tau_el
                Dif = vF2 * tau_el / 3.0
                Lphi = np.sqrt(Dif * tau_phi)
                dQ = np.array([np.linalg.norm(Bmat[i]) / shape[i]
                               for i in range(3)])
                qcart = np.array([np.linalg.norm(q @ Bmat)
                                  for q in Qf if np.linalg.norm(q) > 1e-9])
                Lsc = 2.0 * np.pi / qcart.min() if qcart.size else np.inf
                logger.info("  Cooperon scales: l_el {:.1f} A, L_phi "
                            "{:.1f} A, supercell {:.1f} A, dQ_max {:.3g} "
                            "1/A".format(ell, Lphi, Lsc, dQ.max()))
                if dQ.max() > 1.0 / max(Lphi, 1e-9):
                    need = [max(int(np.ceil(np.linalg.norm(Bmat[i])
                                            * Lphi)), shape[i])
                            for i in range(3)]
                    logger.warning(
                        "  Q-mesh too coarse for the dephasing length: "
                        "dQ_max {:.3g} > 1/L_phi {:.3g} 1/A -- the MESH, "
                        "not the physics, cuts the Cooperon log. Use "
                        "--kmesh >= {} (N_i >= |b_i| * L_phi)."
                        .format(dQ.max(), 1.0 / Lphi, need))
                if ell > Lsc:
                    logger.warning(
                        "  elastic mean free path l_el {:.1f} A exceeds "
                        "the supercell linear size {:.1f} A: periodic "
                        "disorder replicas are phase-coherent and "
                        "contaminate the crossed channel. Minimal "
                        "supercell linear size ~ l_el (increase det S)."
                        .format(ell, Lsc))
                if setup.nc < 8 or asym > 0.05:
                    logger.warning(
                        "  Cooperon with {} configurations / parity "
                        "asymmetry {:.1%}: the Q -> 0 pole region "
                        "amplifies ensemble noise; treat the crossed "
                        "contribution qualitatively and enlarge the "
                        "ensemble.".format(setup.nc, asym))
                kpts_arr = np.asarray(self.kpoints, float)
                n3 = np.asarray(shape, int)

                def _mesh_idx(kf):
                    x = np.mod(kf, 1.0) * n3
                    ix = np.rint(x).astype(int) % n3
                    return int(grid2lrt[ix[0] * n3[1] * n3[2]
                                        + ix[1] * n3[2] + ix[2]])

                def _coop_worker(_st, iQ):
                    Qv = kpts_arr[int(iQ)]
                    part = {}
                    for orb, targ in orb_targets:
                        for eb in targ:
                            res = _cooperon_orbitQ(
                                orb, Qv, kpts_arr, shape, jshift, dVband,
                                wc, a, Gt, Zw, vband, eb * ebin,
                                _mesh_idx)
                            for ik, C in res.items():
                                key = (ik, eb)
                                part[key] = part.get(key, 0.0) + C
                    return part
                tile_map = getattr(ds, 'tile_map', None)
                if tile_map is not None:
                    parts = tile_map(_coop_worker, None,
                                     list(range(self.nkp)))
                else:
                    parts = [_coop_worker(None, iQ)
                             for iQ in range(self.nkp)]
                for part in parts:
                    for key, C in part.items():
                        Cmap[key] = Cmap.get(key, 0.0) + C
                cmag = max((float(np.abs(C).max())
                            for C in Cmap.values()), default=0.0)
                logger.info("  Cooperon vertex accumulated over {} Q "
                            "points; max |C| = {:.3e}"
                            .format(self.nkp, cmag))

            # dressing pass: ladder solve per (orbit, bin) + Cooperon add
            for orb, targ in orb_targets:
                pos = {ik: i for i, ik in enumerate(orb)}
                for eb, items in targ.items():
                    if ladder_on:
                        Lam = _solve_vertex_orbit(
                            orb, jshift, dVband, wc, a, Gt, Zw, vband,
                            eb * ebin, svd_cut, stats)
                    for (ik, n, m) in items:
                        i = pos[ik]
                        v3 = vband[ik][n, m]           # (3,) complex
                        L3 = (Lam[:, i, n, m].copy() if ladder_on
                              else v3.copy())
                        if coop_on and (ik, eb) in Cmap:
                            L3 = L3 + Cmap[(ik, eb)][:, n, m]
                        raw = np.conj(v3[dir1]) * L3[dir2]
                        if self.ortho:
                            M6 = raw[:3].real
                        else:
                            M6 = np.empty(9)
                            M6[:6] = raw.real
                            M6[6:] = raw[:3].imag
                        if n == m:
                            MdiagD[iTl, ispin, ik, n] = M6
                            vertD[iTl, ispin, ik, n] = L3
                        else:
                            inter[(iTl, ispin, ik)][(n, m)] = (L3, M6)
            # ratios (vectorised, relative guard per (band, dir))
            vd = np.stack([vband[ik][np.arange(norb), np.arange(norb)]
                           for ik in range(self.nkp)])   # (nkp,norb,3)
            bare = np.real(np.conj(vd[..., dir1])
                           * vd[..., dir2])[..., :ndir]
            ref_scale = np.abs(bare).max(axis=0, keepdims=True)
            good = np.abs(bare) > 1e-6 * np.maximum(ref_scale, 1e-300)
            lamD[iTl, ispin] = np.where(
                good, MdiagD[iTl, ispin] / np.where(good, bare, 1.0), 1.0)
            guard_skipped += int(np.sum(~good))
            lam3 = lamD[iTl, ispin][..., :3]
            logger.info(
                "  vertex[{} | {}]: lambda(diag) min {:.3f} max {:.3f} "
                "mean {:.3f}; interband pairs {}; dropped modes {} "
                "(cond {:.1e}); source leak {:.1e}{}".format(
                    vertex_channels, vertex_eval, lam3.min(), lam3.max(),
                    lam3.mean(), stats["pairs"], stats["dropped"],
                    stats["cond"], stats["leak"],
                    "  [WARNING: leak > 1e-3, parity protection broken; "
                    "lambda unreliable near affected energies]"
                    if stats["leak"] > 1e-3 else ""))

    if guard_skipped:
        logger.info("  lambda ratio guard: {} (k,band,dir) entries with "
                    "negligible bare moment set to lambda=1 "
                    "(dressed values are still exact)".format(
                        guard_skipped))
    if not per_t and vca:
        pass                                            # shifts already set
    if not per_t:
        for jt in range(1, nt):
            bshift_out[jt] = bshift_out[0]
    elif not combined:
        for jt in range(1, nt):
            bshift_out[jt] = bshift_out[0]

    # ── scattering file ───────────────────────────────────────────────────
    gamma_out = np.maximum(gamma_out, 1e-14)
    scatobj.defineScatteringRates(gamma_out, qpweight=Z_out,
                                  bandshift=bshift_out)
    mt0 = os.path.getmtime(scat_output) if os.path.isfile(scat_output) \
        else None
    scatobj.createOutput(scat_output)
    if os.path.isfile(scat_output) and (mt0 is None or
                                        os.path.getmtime(scat_output) > mt0):
        logger.info("computeDisorder: scattering -> {}".format(scat_output))

    # ── dressed energy file(s) ────────────────────────────────────────────
    if vertex and dressed_output:
        t_write = list(range(nt)) if (dressed_per_t and per_t) else [0]
        base, ext = os.path.splitext(dressed_output)
        for jT in t_write:
            fout = dressed_output if len(t_write) == 1 else \
                '{}_T{:04d}{}'.format(base, jT + 1, ext)
            shutil.copyfile(energy_hdf5, fout)
            with h5py.File(fout, 'r+') as h:
                h.attrs['disorder_dressed'] = str(
                    'channels=' + vertex_channels + ', vertex_eval='
                    + vertex_eval + ', rate_eval=' + rate_eval)
                h.attrs['disorder_archive'] = str(archive_path)
                h.attrs['temperature_K'] = float(Ts[jT])
                if gamma_add != (0.0, 0.0):
                    h.attrs['gamma_add_c0_c2'] = np.asarray(gamma_add)
                for isp in range(spins):
                    pre = '' if spins == 1 else ('up/' if isp == 0
                                                 else 'dn/')
                    h[pre + 'momentsDiagonal'][...] = MdiagD[jT, isp]
                    for (kT, jsp, ik), dd in inter.items():
                        if kT != jT or jsp != isp:
                            continue
                        key = pre + 'kPoint/{:010d}/moments'.format(ik + 1)
                        if key not in h:
                            continue
                        M = h[key][()]
                        for (n, m), (_L3, M6) in dd.items():
                            M[n, m, :] = M6
                            M[m, n, :] = M6
                        for n in range(norb):
                            M[n, n, :] = MdiagD[jT, isp, ik, n]
                        h[key][...] = M
            logger.info("computeDisorder: dressed energy -> {}".format(
                fout))
        if per_t and nt > 1 and not dressed_per_t:
            logger.warning(
                "vertex dressing is temperature dependent (combined mode "
                "or --gamma-add) but only the FIRST temperature's dressed "
                "file was written; use --dressed-per-t for one file per "
                "temperature (then run linretrace per T).")

    # ── lambda analysis file ──────────────────────────────────────────────
    if vertex and lambda_output:
        with h5py.File(lambda_output, 'w') as h:
            h.attrs['identifier'] = str('LRTC-disorder-lambda')
            h.attrs['archive'] = str(archive_path)
            h.attrs['vertex_eval'] = str(vertex_eval)
            h.attrs['rate_eval'] = str(rate_eval)
            h.attrs['vertex_channels'] = str(vertex_channels)
            h.attrs['gamma_add_c0_c2'] = np.asarray(gamma_add)
            h['.kmesh/points'] = np.asarray(self.kpoints, float)
            h['.kmesh/nkp'] = self.nkp
            h['.quantities/mu'] = mu
            h['.quantities/tempAxis'] = Ts
            for isp in range(spins):
                pre = '' if spins == 1 else ('up/' if isp == 0 else 'dn/')
                h[pre + 'lambdaDiagonal'] = (lamD[:, isp] if per_t
                                             else lamD[0, isp])
                h[pre + 'vertexDiagonal'] = (vertD[:, isp] if per_t
                                             else vertD[0, isp])
                h[pre + 'momentsDiagonalDressed'] = (
                    MdiagD[:, isp] if per_t else MdiagD[0, isp])
                for (jT, jsp, ik), dd in inter.items():
                    if jsp != isp or (not per_t and jT != 0):
                        continue
                    arr = np.zeros((len(dd), 8))
                    for i, ((n, m), (L3, _M6)) in enumerate(dd.items()):
                        arr[i, 0], arr[i, 1] = n + 1, m + 1
                        arr[i, 2:5] = np.real(L3)
                        arr[i, 5:8] = np.imag(L3)
                    key = (pre + ('T{:04d}/'.format(jT) if per_t else '')
                           + 'kPoint/{:010d}/lambdaPairs'.format(ik + 1))
                    h[key] = arr
            h['correlatorW2'] = V2_export
            h['correlatorW2'].attrs['note'] = (
                'ensemble-averaged squared correlator '
                'sum_c w_c |dV^band_c|^2, shape (nkp, N_transfers, nb, '
                'nb): element [k, p, n, m] couples state (k, n) to '
                '(k - q_p, m); band basis of the LAST processed '
                '(spin, T) pass')
            h['bandEnergies'] = a_export
            h['railWidths'] = Gt_export
            h.attrs['note'] = str(
                'lambdaDiagonal = dressed/bare diagonal transition matrix '
                'elements; vertexDiagonal = Lambda^alpha_nn (complex); '
                'lambdaPairs rows: n, m, Re Lambda^{x,y,z}, Im '
                'Lambda^{x,y,z}. Guarded entries (negligible bare moment) '
                'carry lambda=1.')
        logger.info("computeDisorder: lambdas -> {}".format(lambda_output))

    logger.info("computeDisorder: done.")
    return dict(gamma=gamma_out, Z=Z_out, bandshift=bshift_out,
                lambda_diag=lamD, moments_dressed=MdiagD)
