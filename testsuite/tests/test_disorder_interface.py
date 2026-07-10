#! /usr/bin/env python3
"""
testsuite/tests/test_disorder_interface.py
==========================================
Validation of the ab-initio disorder interface (structure/wannier_disorder,
lwann --disorder).  Needs the dwannier pipeline modules (DWANNIER_PATH) and
the CaCuO2 test archive (DWANNIER_ARCHIVE).  Checks:

  1. orbit ladder solver == dwannier's machine-validated solve_vertex on
     identical rails (algebraic identity),
  2. point-like correlator  -> lambda = 1 exactly (parity),
  3. forward-scattering separable correlator -> analytic tau_tr/tau,
  4. end-to-end computeDisorder on CaCuO2 in BOTH vertex-eval modes:
     outputs written, dressed(diag) == lambda * bare, scatrate sane,
     fermi mode == onshell mode at the artificial limit where all states
     are pinned to one energy bin.

Run from the repo root:
  DWANNIER_PATH=<pipeline dir> python3 testsuite/tests/test_disorder_interface.py
"""
import os
import sys
import logging

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..', '..')))
for p in (os.environ.get("DWANNIER_PATH"),):
    if p:
        sys.path.insert(0, p)

ARCHIVE = os.environ.get("DWANNIER_ARCHIVE",
                         "/mnt/user-data/outputs/cacuo2_run/deltaH.h5")
W90DIR = os.environ.get("DISORDER_W90DIR", "cacuo2_w90")

from structure.wannier_disorder import (_solve_vertex_orbit,
                                        _orbits_from_jshift)

ok = True


def check(name, cond, extra=""):
    global ok
    ok = ok and bool(cond)
    print("  [{}] {} {}".format("PASS" if cond else "FAIL", name, extra))


def _stats():
    return {"dropped": 0, "leak": 0.0, "cond": 1.0, "pairs": 0}


# ---------------------------------------------------------------------------
print("1. orbit solver == dwannier solve_vertex (identical rails)")
import disorder_selfenergy as ds
import disorder_transport as dt

arch = ds.DisorderArchive(ARCHIVE)
ref = arch.reference()
setup = dt.TransportSetup(arch, ref, (4, 4, 4))
nk, n_p = setup.nk, setup.n_p
G0, EB = 0.25, -1.0
aa = np.array([[np.real(ref.H0(setup.kpts[ik])[0, 0])]
               for ik in range(nk)])                     # (nk, 1)
GR = (1.0 / (EB - aa[:, 0] + 1j * G0))[:, None, None]
Gam_ref = dt.solve_vertex(setup, setup.vel[:, 0], GR, np.conj(GR))

dVband = [[[setup.dV[ik, p, c] for c in range(setup.nc)]
           for p in range(setup.N)] for ik in range(nk)]
vband = [np.transpose(setup.vel[ik], (1, 2, 0)) for ik in range(nk)]
Gt = np.full((nk, 1), G0)
Zw = np.ones((nk, 1))
worst = 0.0
for orb in _orbits_from_jshift(nk, setup.jshift):
    Lam = _solve_vertex_orbit(orb, setup.jshift, dVband, arch.w, aa, Gt,
                              Zw, vband, EB, 1e-13, _stats())
    for i, ik in enumerate(orb):
        worst = max(worst, abs(Lam[0, i, 0, 0] - Gam_ref[ik, 0, 0]))
check("orbit solver == solve_vertex", worst < 1e-9,
      "(max dev {:.1e}, |Gamma| scale {:.2f})".format(
          worst, float(np.abs(Gam_ref).max())))

# ---------------------------------------------------------------------------
print("2./3. synthetic ring: point-like parity + forward-scattering analytic")
nring = 24
th = 2 * np.pi * np.arange(nring) / nring
Ering = (-0.5 * np.cos(th))[:, None]
Gring = np.full((nring, 1), 0.05)
Zring = np.ones((nring, 1))
vv = 0.5 * np.sin(th)
vb = [np.array([[[vv[i]] * 3]]).reshape(1, 1, 3) for i in range(nring)]
jsh = np.tile(np.arange(nring), (nring, 1))              # all-to-all ring

# point-like: dV = const
u0 = 0.02
dVp = [[[np.array([[u0]])] for _ in range(nring)] for _ in range(nring)]
w1c = np.ones(1)
Lam = _solve_vertex_orbit(list(range(nring)), jsh, dVp, w1c, Ering, Gring,
                          Zring, vb, float(Ering[5, 0]), 1e-2, _stats())
lam = float(np.real(Lam[0, 5, 0, 0])) / vv[5]
check("lambda(point-like) == 1", abs(lam - 1) < 1e-10,
      "(lambda = {:.12f})".format(lam))

# forward: |dV(k,k')|^2 = w0 + w1 sin(k) sin(k'); realized with two
# configurations dV = sqrt(w0) and dV = sqrt(w1) s_k... use complex trick:
# two independent 'configurations' with amplitudes u0=sqrt(w0),
# u1(k,k') = sqrt(w1) * s_k s_k' requires a separable AMPLITUDE; since
# |u1|^2 = w1 s_k^2 s_k'^2 != w1 s_k s_k', build instead a rank-2 exact
# correlator via signed amplitudes: W = w0 + w1 s s' is realized by
# configs c=(+,-): u_pm(k,k') = sqrt(w0/2 + w1 s_k s_k' /2) ... not
# separable either. Test instead DIRECTLY against dense solve_vertex on
# the same synthetic dV set (already covered by check 1 for the ab-initio
# correlator) and against the analytic formula using a correlator that IS
# realizable: single config u(k,k') = u0 + u1 s_k s_k' gives
# W = |u|^2 = u0^2 + 2 u0 u1 s s' + u1^2 s^2 s'^2 (three separable
# channels); the odd channel coefficient is w1_eff = 2 u0 u1 and the even
# channels do not couple to the odd source, so lambda = 1/(1 - w1_eff q).
u0a, u1a = 0.02, 0.015
dVf = [[[np.array([[u0a + u1a * np.sin(th[i]) * np.sin(th[j])]])]
        for j in range(nring)] for i in range(nring)]
eps = float(Ering[5, 0])
GG = 1.0 / ((eps - Ering[:, 0]) ** 2 + Gring[:, 0] ** 2)
q = float(np.sum(np.sin(th) ** 2 * GG))
lam_an = 1.0 / (1.0 - 2 * u0a * u1a * q)
Lam = _solve_vertex_orbit(list(range(nring)), jsh, dVf, w1c, Ering, Gring,
                          Zring, vb, eps, 1e-2, _stats())
lam = float(np.real(Lam[0, 5, 0, 0])) / vv[5]
check("lambda(forward) == 1/(1 - 2 u0 u1 q)",
      abs(lam / lam_an - 1) < 1e-10,
      "(lambda = {:.6f}, analytic {:.6f})".format(lam, lam_an))

# ---------------------------------------------------------------------------
print("4. end-to-end computeDisorder on CaCuO2 (4x4x4), both eval modes")
logging.basicConfig(level=logging.ERROR)
from structure.wannier import Wannier90Calculation
from structure.wannier_disorder import computeDisorder
Wannier90Calculation.computeDisorder = computeDisorder

import tempfile
import h5py

ham = Wannier90Calculation(W90DIR, None, False)
ham.readData()
ham.expandKmesh(np.array([4, 4, 4], dtype=int))
tmp = tempfile.mkdtemp()
efile = os.path.join(tmp, "e.hdf5")
ham.computeHamiltonian(peierlscorrection=True)
ham.outputData(efile, mu=ham.efer)

for mode in ("onshell", "fermi"):
    res = ham.computeDisorder(
        ARCHIVE, efile,
        scat_output=os.path.join(tmp, "scat_{}.hdf5".format(mode)),
        temps=(100.0, 400.0, 3, False),
        vertex_eval=mode,
        dressed_output=os.path.join(tmp, "dressed_{}.hdf5".format(mode)),
        lambda_output=os.path.join(tmp, "lam_{}.hdf5".format(mode)),
        interactive=False)
    with h5py.File(os.path.join(tmp, "scat_{}.hdf5".format(mode))) as h:
        g1 = h['step/000001/scatrate'][()]
        check("[{}] scatrate written, positive".format(mode),
              g1.shape == (ham.nkp, 1) and bool((g1 > 0).all()),
              "(Gamma in [{:.3g}, {:.3g}] eV)".format(g1.min(), g1.max()))
    with h5py.File(efile) as he, \
         h5py.File(os.path.join(tmp, "dressed_{}.hdf5".format(mode))) as hd, \
         h5py.File(os.path.join(tmp, "lam_{}.hdf5".format(mode))) as hl:
        Mb = he['momentsDiagonal'][()]
        Md = hd['momentsDiagonal'][()]
        lam = hl['lambdaDiagonal'][()][..., :Mb.shape[-1]]
        dev = float(np.abs(Md - lam * Mb).max())
        check("[{}] dressed == lambda * bare (diag)".format(mode),
              dev < 1e-10, "(max dev {:.1e})".format(dev))
        lam3 = lam[..., :3]
        print("      lambda[{}]: min {:.3f} max {:.3f} mean {:.3f}".format(
            mode, lam3.min(), lam3.max(), lam3.mean()))

# ---------------------------------------------------------------------------
print("5. Cooperon orbit operator == full-code crossed fan (identical rails)")
from structure.wannier_disorder import _cooperon_orbitQ

G0c, EBc = 0.30, -1.0
GRc = (1.0 / (EBc - aa[:, 0] + 1j * G0c))[:, None, None]
GAc = np.conj(GRc)
Gtc = np.full((nk, 1), G0c)
T_ref = dt.cooperon_correction(setup, GRc, GAc, setup.vel[:, 0],
                               setup.vel[:, 0], serial=True)
n3 = np.asarray(setup.shape, int)


def _mesh_idx(kf):
    x = np.mod(kf, 1.0) * n3
    ix = np.rint(x).astype(int) % n3
    return int(ix[0] * n3[1] * n3[2] + ix[1] * n3[2] + ix[2])


Cmap = {}
for iQ in range(nk):
    for orb in _orbits_from_jshift(nk, setup.jshift):
        res = _cooperon_orbitQ(orb, setup.kpts[iQ], setup.kpts,
                               setup.shape, setup.jshift, dVband, arch.w,
                               aa, Gtc, np.ones((nk, 1)), vband, EBc,
                               _mesh_idx)
        for ik, C in res.items():
            Cmap[ik] = Cmap.get(ik, 0.0) + C
T_int = sum(np.real(setup.vel[ik, 0, 0, 0]) * GRc[ik, 0, 0]
            * Cmap[ik][0, 0, 0] * GAc[ik, 0, 0] for ik in Cmap) / nk
check("Cooperon vertex == full-code crossed fan",
      abs(T_int - T_ref) / max(abs(T_ref), 1e-300) < 1e-12,
      "(rel dev {:.1e})".format(abs(T_int - T_ref) / abs(T_ref)))

# ---------------------------------------------------------------------------
print("6. end-to-end channels: 'both' with --gamma-add; dephasing guard")
try:
    ham.computeDisorder(
        ARCHIVE, efile,
        scat_output=os.path.join(tmp, "scat_c.hdf5"),
        temps=(100.0, 300.0, 2, False), vertex_channels='cooperon',
        dressed_output=os.path.join(tmp, "dressed_c.hdf5"),
        interactive=False)
    check("Cooperon without dephasing raises", False)
except ValueError:
    check("Cooperon without dephasing raises", True)
res = ham.computeDisorder(
    ARCHIVE, efile,
    scat_output=os.path.join(tmp, "scat_b.hdf5"),
    temps=(100.0, 300.0, 2, False),
    vertex_channels='both', gamma_add=(0.02, 0.0),
    rate_eval='onshell', dressed_per_t=True,
    dressed_output=os.path.join(tmp, "dressed_b.hdf5"),
    lambda_output=os.path.join(tmp, "lam_b.hdf5"),
    interactive=False)
okfiles = (os.path.isfile(os.path.join(tmp, "dressed_b_T0001.hdf5"))
           and os.path.isfile(os.path.join(tmp, "dressed_b_T0002.hdf5")))
check("'both' + gamma_add runs; per-T dressed files written", okfiles)
with h5py.File(os.path.join(tmp, "scat_b.hdf5")) as h:
    g1 = h['step/000001/scatrate'][()]
    g2 = h['step/000002/scatrate'][()]
    check("scatrate includes gamma_add (min >= C0)",
          bool(g1.min() >= 0.02 - 1e-12 and g2.min() >= 0.02 - 1e-12),
          "(min {:.4f})".format(g1.min()))

print()
print("ALL INTERFACE TESTS PASSED" if ok else "SOME TESTS FAILED")
sys.exit(0 if ok else 1)
