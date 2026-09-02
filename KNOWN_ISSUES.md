# Known issues, observations and deferred work

Register of defects, latent hazards and performance opportunities that have
been *identified* but not (yet) resolved. The purpose is to make sure nothing
found during development is lost simply because attention moved elsewhere.

Each entry is written to be actionable on its own: what the symptom is, where
in the tree it lives, why it happens, how it was observed, and what a fix would
have to do. Someone picking up an entry cold should not need to re-derive the
analysis.

**Conventions**

- Entries are numbered `KI-nn` (issues) and `PERF-nn` (performance work) and are
  never renumbered or deleted; resolved entries move to
  [Resolved](#resolved) with the commit or patch that closed them.
- *Severity*: **high** = wrong physical results possible; **medium** = wrong
  results only in a restricted regime, or a hazard that is currently latent;
  **low** = cosmetic, ergonomic, or hygiene.
- *Confidence*: whether the diagnosis is confirmed by a reproducer or is still
  a reading of the code.

Last updated: 2026-09-02 (patch 20-21).

---

## Open issues

### KI-01 — Eigenvector phase makes off-diagonal `momentsBfield` gauge dependent

| | |
|---|---|
| Severity | **high** (raised 2026-08-08) |
| Confidence | confirmed (two reproducers below); consumption confirmed in `response.F90` |
| Location | `structure/wannier.py:computeHamiltonian`, `structure/tb.py:_computeHk` |
| Affects | `/kPoint/*/momentsBfield` (inter-band only) -> `L11B`/`L12B`/`L22B` inter-band contributions, i.e. runs with `BFieldMode` **and** full B moments |
| Not affected | energies, `moments`, `momentsDiagonalBfield`, and every band-diagonal magnetic quantity |

**Symptom.** Band off-diagonal B-field moments change when the k-point batch
size changes, even though the physics does not. Diagonal elements, energies and
`moments` are unaffected.

**Cause.** `np.linalg.eigh` fixes eigenvectors only up to a phase
`U_n -> e^{i phi_n} U_n` (a sign, for real symmetric input). LAPACK does not
guarantee a particular choice, and the choice can depend on the batch size or
the LAPACK build. Matrix elements in the band basis inherit the phase:

- `moments`: `(v^i_nm)^* v^j_nm` — the phases cancel, gauge invariant.
- `momentsBfield`: `eps_cij (v^a_nm)^* v^i_nm c^bj_nm` — three factors carrying
  the *same* pair `(n,m)`, so one factor of `e^{i(phi_m - phi_n)}` survives.
  Gauge dependent for `n != m`; invariant for `n == m`.

This is pre-existing and independent of the k-blocking introduced in patch 05 —
blocking merely made it visible, because it changes the batch handed to LAPACK.

**Escalation: on the symmetrising path the magnitude is arbitrary too.**
The optical / B-field moments of an irreducible k-point are the average over
its star, `M_sym = (1/nsym) sum_S M(S k)`. Each term carries its *own*,
independently arbitrary phase `e^{i(phi_m(Sk) - phi_n(Sk))}`. Summing complex
numbers with mutually random phases does not merely rotate the result — it
changes its modulus. So for `n != m` the star-averaged `momentsBfield` is not
gauge covariant, it is simply not well defined.

Measured by running the *same* input through the pre- and post-patch-08
symmetrising code (which differ only in the batch shape handed to `eigh`, and
therefore only in the phases LAPACK returns), 584 k-points, `nsym = 48`:

```
|Bopt| off-diagonal, scale        : 9.35e-03
|Bopt| off-diagonal, deviation    : 2.73e-03   (29 % relative)
|Bopt| band-diagonal, deviation   : 8.33e-17   (exact)
```

418 of 584 k-points are affected, with band gaps from 0.16 to 3.0 eV — this is
not a degenerate-subspace effect. The band-diagonal elements, `moments`, and
the energies are all unaffected, so this is confined to
`momentsBfield[n != m]`.

Note the consequence for reproducibility: those numbers can change with the
LAPACK build, the NumPy version, or the block size, with no warning.

**Reproducer.**

```python
h_a = run(kblock=None)   # single pass over 64 k-points
h_b = run(kblock=1)
np.abs(np.abs(h_a.velocities[0]) - np.abs(h_b.velocities[0])).max()   # 2.2e-16
h_b.velocities[0][6,0,1,:] / h_a.velocities[0][6,0,1,:]
#   -> -1.0 + 1.1e-17j   (|ratio| == 1: a pure phase, here a sign flip)
```

The band gap at that k-point is 0.8 eV, so this is *not* a degenerate-subspace
ambiguity — it is the ordinary U(1) freedom of a non-degenerate eigenvector.

**Is it used? Yes, and linearly.** (`src_linretrace/response.F90`, read
2026-08-08; nothing in the Fortran core was changed.)

The inter-band magnetic response is accumulated as

```fortran
if (algo%lBfield .and. edisp%lBFullMoments) then
  do idir1=1,3 ; do idir2=1,3 ; do idir3=1,3
    resp%sB_full(idir3,idir2,idir1,iband1,is,info%ik) = &
    resp%sB_full(idir3,idir2,idir1,iband1,is,info%ik) + calc_sigmaB &
                 * edisp%MBopt(idir3,idir2,idir1,iband2,iband1,is,info%ik)
```

inside the `iband1`/`iband2` double loop, i.e.

    sB_full[..., n] = sum_m  K(eps_n, eps_m, Gamma_n, Gamma_m) * MBopt[..., m, n]

`calc_sigmaB` (and `calc_alphaB`, `calc_xiB`) are built only from band
energies, quasi-particle weights and scattering rates -- all gauge invariant --
and enter as a **linear** coefficient. So each term in the sum over `m` carries
its own arbitrary factor `e^{i(phi_m - phi_n)}`, the phases do not cancel, and
the sum is arbitrary. `L11B`, `L12B`, `L22B` inter-band contributions inherit
this directly.

Contrast the non-magnetic response a few lines above, which uses
`edisp%Mopt(idir,iband2,iband1)`. `Mopt` is `v_a^*[n,m] v_b[n,m]` -- bilinear
in the *same* ordered pair -- so the phases cancel inside each element and
`resp%s_full` is gauge invariant. The band-diagonal magnetic path
(`MBoptDiag`, `n == m`) is likewise invariant and is the one that is
unconditionally required (`if (algo%lBfield .and. .not. edisp%lBIntraMoments)
call stop_with_message`). Only the optional `lBFullMoments` path is affected.

**Structural observation (hypothesis, not a claim).** The moment is built as

    MBopt[n,m] ~ eps_zij  v_a^*[n,m]  v_i[n,m]  c_yj[n,m]

-- three factors all carrying the *same* ordered band pair `(n,m)`, hence a net
`e^{i(phi_m - phi_n)}`. A quantity built from band-basis matrix elements is
gauge invariant iff its band indices form closed loops: `v^*[n,m] v[m,n]`
(two-index loop, as in `Mopt`), or `v[n,m] v[m,l] c[l,n]` (three-index loop).
The present expression is not of either form, which suggests the trilinear
object should contract over an intermediate band index rather than repeat the
pair. That is consistent with the expectation that the correct multi-orbital
expression is gauge invariant -- and with the possibility that terms are
missing: for `n != m` the plain curvature matrix element `c = d2H/dk dk`
rotated into the band basis is not the gauge-covariant derivative of the
inter-band velocity, which additionally involves the Berry connection
`A[n,m] = i <u_n | d_k u_m>`. This is exactly the kind of term a careful
derivation of the response at finite vector potential and finite `q` with
orbital indices would settle.

**A cheap consistency test for any candidate expression.** Apply a random
diagonal phase `U -> U diag(e^{i phi_n})` after `eigh` and recompute. The
correct expression must be invariant to machine precision. The current one is
not: see the 29 % magnitude shift above, which arose merely from changing the
batch shape handed to LAPACK.

**Open questions before fixing.**

1. Settle the derivation first. Gauge-fixing the eigenvectors would make the
   present numbers *reproducible* without making them *right*, and would
   therefore mask the problem rather than solve it. If a gauge convention is
   introduced anyway -- and it is worth having for bitwise reproducibility
   across LAPACK/NumPy versions -- it should be introduced together with a
   note that it does not validate `lBFullMoments` output.
2. If a fix is needed, the standard remedy is a gauge fixing convention applied
   right after `eigh`, e.g. rotate each eigenvector so that its
   largest-modulus component is real and positive. That is cheap
   (`O(nkp * nproj^2)`), deterministic, and would also make output bitwise
   reproducible across NumPy/LAPACK versions — a desirable property in its own
   right. It must be applied identically in `tb.py` and `wannier.py`.
3. Note that a *global* gauge fixing does not make inter-band quantities
   physically meaningful if they are not gauge invariant to begin with; it only
   makes them reproducible.

---

### KI-02 — `momentsBfield` is written per k-point as an individual HDF5 dataset

| | |
|---|---|
| Severity | low (correctness) / high (practical) |
| Confidence | confirmed |
| Location | `structure/inout.py:h5output` |
| Blocked by | the Fortran reader (`src_linretrace/input.f90`) |

**Symptom.** Writing dominates the runtime of a refined-mesh generator run. On
a 2-orbital test at 5056 k-points: compute 0.23 s, output **3.23 s**. The cost
scales with the number of k-points and is essentially independent of `nproj`,
i.e. it is per-group overhead rather than data volume.

**Cause.** When `opticfull` / `bopticfull` are set, `h5output` creates one
dataset per k-point:

```python
for ikp in range(escalc.nkp):
    h5out[prefix+'kPoint/{:010}/moments'.format(ikp+1)] = ...
```

For an adaptively refined mesh with 1e5–1e6 k-points this is 1e5–1e6 group and
dataset creations.

**Why it is not fixed.** The obvious remedy — a single chunked
`(nkp, nbands, nbands, ndir)` dataset — changes the file layout and would
require a matching change in the Fortran reader. The Fortran core is currently
held fixed by policy.

**Workaround in place.** `--intraonly` skips the per-k-point datasets entirely
(output 3.23 s -> 0.02 s on the same case), at the cost of discarding
inter-band elements.

**If it is ever revisited.** A backwards-compatible route exists: write the
batched dataset *in addition*, mark it with an attribute, and have the reader
prefer it when present. Old files keep working, new files are fast.

---

### KI-03 — `structure/wannier.py:readData` refuses irreducible Wannier90 grids

| | |
|---|---|
| Severity | low |
| Confidence | confirmed (explicit `raise`) |
| Location | `structure/wannier.py:readData` |

`readData` raises `IOError('Irreducible momentum grids not implemented')`
whenever `nkx*nky*nkz != nkp`. Only the `--wien2k` path produces an irreducible
Wannier mesh.

Consequence for the refinement: `lwann refine` gains the symmetrisation
machinery of patch 06, but in the normal path there is nothing to symmetrise,
because the coarse mesh is always the full BZ. The k-point saving from
symmetry — the main reason the machinery exists — is therefore only available
to `ltb` and to `--wien2k`-derived Wannier runs.

Lifting this would mean deriving the symmetry operations for a Wannier90 mesh
(spglib on the `.win` unit cell, as `ltb` does) and folding the grid, which is
a self-contained piece of work.

---

### KI-04 — A custom (refined) mesh silently breaks any grid-indexed code path

| | |
|---|---|
| Severity | medium (latent) |
| Confidence | confirmed by inspection |
| Location | `structure/wannier_disorder.py` (out-of-tree), any future grid-indexed consumer |

`setCustomKmesh` sets `nkx = nky = nkz = 1` by convention, because a refined
mesh has no regular grid spacing. Code that reconstructs a grid address from
k-points will then compute nonsense rather than fail:

```python
shape = (self.nkx, self.nky, self.nkz)          # (1,1,1) for a custom mesh
lin   = _mesh_index_map(self.kpoints, shape)    # every k-point maps to index 0
```

The disorder extension does exactly this (`computeDisorder`, and the Cooperon
`_mesh_idx` helper). It is not distributed with the public branch, so nothing in
the repository currently trips over it, but the hazard is real for anyone
combining the two.

**Hook already in place.** `Wannier90Calculation.customMesh` (and
`TightBinding.customMesh`) is set to `True` by `setCustomKmesh`. Any grid-indexed
routine should begin with

```python
if getattr(self, 'customMesh', False):
    raise NotImplementedError('<routine> requires a regular k-grid; '
                              'this object carries a custom/refined mesh.')
```

---

### KI-05 — `structure/units.py` is not the single source of truth for constants

| | |
|---|---|
| Severity | low |
| Confidence | confirmed |
| Location | see table |

The Boltzmann constant is defined independently in five places, with three
different values:

| location | value |
|---|---|
| `src_linretrace/params.f90` | `8.6173324e-5` |
| `scattering/dmftrenormalization.py` | `8.61733034e-5` |
| `scattering/fullscattering.py` | `8.61733034e-5` |
| `scripts/selfconsistent_chemicalpotential.py` | `8.61733034e-5` |
| `postproc/output.py` (`_kB_eV`) | `8.617333262e-5` |

Patch 01 added `kB_eV = 8.61733034e-5` to `structure/units.py` and used it in
`structure/meshrefine.py`, but did not touch the existing local definitions
(they are in working code). The differences are ~1e-8 relative and therefore
harmless numerically; the hazard is that a *sixth* copy gets added, or that one
is edited and the others are not. The origin of this class of bug was a missing
constant, not a wrong one — see [Resolved](#resolved), KI-R1.

KI-R9 is the same failure mode on the Fortran side, one severity level up: a
constant that was typed rather than derived, and typed wrongly. Both argue for
deriving constants from a single definition wherever the arithmetic allows it.

Suggested cleanup: import `kB_eV` from `structure.units` everywhere on the
Python side and adopt the CODATA-2018 value `8.617333262e-5` consistently. The
Fortran parameter should be aligned in the same pass, but that is a Fortran-core
change.

---

### KI-06 — `momentsBfield` and `velocities` are not populated on the irreducible path

| | |
|---|---|
| Severity | low |
| Confidence | confirmed by inspection |
| Location | `structure/wannier.py:computeHamiltonian`, irreducible branch |

The irreducible branch appends `energies`, `opticalMoments` and
`BopticalMoments` but never `velocities` / `curvatures`; the reducible branch
appends all five. Nothing in `h5output` needs the raw velocities, so this is
currently invisible, but it makes the two branches non-interchangeable for any
consumer that reaches for `calc.velocities` (e.g. a future post-processing
routine, or `structure/wannier_dmft.py` if it is ever extended). Worth either
populating them for symmetry or documenting the asymmetry at the call site.

---

### KI-07 — Coarse-mesh dimensionality cannot be inferred, only inherited

| | |
|---|---|
| Severity | low (mitigated) |
| Confidence | confirmed |
| Location | `structure/es.py:_defineDimensionsFromKpoints` |

`ndim` / `dims` are properties of the *physical setup*, not of the k-mesh. A
layered system can carry a fully three-dimensional k-mesh while the Wannier
projection is effectively two-dimensional (e.g. a single planar `d_x2-y2`
orbital in a cuprate). Inferring the dimensions from the spread of the k-points
would report 3 in that case.

Patch 04 made the generators inherit `.unitcell/dims` from the coarse HDF5, so
the refined file always reproduces whatever the parent run decided. The
k-point inference survives only as a fallback for files that predate the
`dims` dataset, and is documented as unable to handle this case. The entry is
kept open as a reminder that the fallback is a trap if it is ever promoted back
to the default.

---

### KI-08 — Reflection branch is three orders less accurate than the direct branch

| | |
|---|---|
| Severity | low (unreachable in practice) |
| Confidence | confirmed |
| Location | `src_linretrace/digamma/psi_fast.F90`, and the same in the CERNLIB originals |

Measured worst relative error over K = 0..4 against `mpmath` references, in
real quad arithmetic:

```
Re z > 0 (every LinReTraCe argument) : 9.81e-34   (~5 x quad eps)
Re z < 0 (reflection branch)         : 1.13e-31   at K = 4, z = (-2.3, 1.1)
```

The loss comes from the reflection formula itself, where `cot(pi u)` cancels;
it is not a property of the recurrence or the series. The old routine has the
same behaviour and is worse (3.8e-30 at the same point), so `psi_fast` is a
~90x improvement there rather than a regression.

Unreachable in LinReTraCe: every argument is
`z = 1/2 + beta/2pi (Gamma + i Z (E - mu))`, so `Re z >= 1/2 > 0` always. The
branch exists only so the module is a drop-in for the general routine.

Recorded so that nobody re-measures it and concludes the module is broken. If
the branch is ever needed, the fix is to evaluate `cot(pi u)` in a form that
does not cancel near the poles, not to touch the series.

---

### KI-09 — The CERNLIB digamma archive is now dead weight

| | |
|---|---|
| Severity | low |
| Confidence | confirmed |
| Location | `src_linretrace/digamma/`, `src_linretrace/Makefile` |

After patch 11 there are no live calls to `wpsipg` or `wpsipghp`: all 23
remaining references are commented out under `! --- deprecated CERNLIB path`
banners, kept so the old behaviour can be restored by uncommenting. The
archive `cern_digamma.a` is still built and linked for that reason.

Retirement checklist, for when the new path has been exercised on production
runs:

1. Delete the commented blocks in `fermi.f90`, `response.F90`, `root.F90`
   (7 call sites; `git log` will preserve them).
2. Drop `wpsipg.o` and `wpsipg_hp.o` from `cern_digamma.a`. Check first
   whether `cpsipg.f` -- the CERNLIB single-precision wrapper, which calls
   `WPSIPG` -- is referenced by anything; if not, the whole archive and the
   `abend/lenocc/mtlprt/mtlset` support files go with it.
3. Remove `DIGAMMA` from `src_linretrace/Makefile` and the `$(DIGAMMA)`
   prerequisite from the `all` target.
4. `kmax_double` / `kmax_quad` in `digamma/Makefile` disappear with it. The
   new module has no preprocessor knobs at all -- the expansion order is a
   run-time decision (see PERF-R2).

Note that `psi_fast.o` is built by the *parent* Makefile with `$(FC)`, not by
the archive: it is a Fortran module, and `.mod` files are not portable between
compilers, while `digamma/Makefile` builds with `$(FCDG)`.

---

### KI-10 — `psi_fast` duplicates its implementation across kinds

| | |
|---|---|
| Severity | low (hygiene) |
| Confidence | confirmed |
| Location | `src_linretrace/digamma/psi_fast.F90` |

`psi_range_qp` / `psi_range_dp`, `psi0_imag_qp` / `psi0_imag_dp` and their
wrappers are line-for-line the same algorithm with a different kind
parameter, which is the same double/quad duplication the old `.F77` pair had
-- moved, not solved. It is at least now confined to one file, and the
coefficient table is shared (`CD = real(CQ, dp)`), so the two kinds cannot
drift apart numerically.

Unifying would mean either a preprocessor include of the body (back to a
`cpp` pass, which the module deliberately avoids) or parameterised derived
types, which gfortran supports only partially. Worth revisiting when the rest
of the double/quad split in the core is addressed; not worth doing on its own.

---

### KI-11 — Exactly degenerate bands give gauge-dependent `momentsDiagonal`, and a spurious gap

| | |
|---|---|
| Severity | medium (restricted regime: exact band degeneracies on the mesh) |
| Confidence | confirmed by reproducer |
| Location | `structure/tb.py:_computeHk`, `structure/wannier.py:computeHamiltonian` |

Sibling of KI-01, on the band-*diagonal* elements this time, and it does not
need the B-field path. Where two bands are exactly degenerate at a mesh point,
`eigh` returns an arbitrary basis inside the degenerate subspace, and
`|v_nn|^2` is not invariant under that rotation.

Graphene at `K = (2/3, 1/3, 0)` on a mesh whose divisions are a multiple of 3:

```
ltb   (H(K) truncated to exactly 0 by the |h|<1e-14 clamp) : momentsDiagonal = 0
lwann (same model round-tripped through hr.dat, H(K) ~ 1e-8): momentsDiagonal = (0.1875, 0.5625, 0)
```

The 1e-8 asymmetry in `hr.dat` lifts the degeneracy numerically and rotates the
eigenvectors into the physically meaningful basis, which is arguably the
*better* answer — the `ltb` value of zero is an artifact of the clamp. The same
1e-8 also opens a fake `7.3e-08 eV` gap, so `.bands/bandgap/gapped` flips to
`True` and `linretrace` switches root method. On a 12x12x1 coarse graphene mesh
the two files give `L11_xx` of 2.1e+02 and 1.4e+07 respectively.

Refinement masks the symptom, because subdivided points no longer sit exactly on
`K`; it does not fix it. A proper fix needs degenerate-subspace handling
(diagonalise the velocity operator within each degenerate block) rather than
relying on the numerical noise to do it. Until then: prefer meshes that avoid
placing points exactly on a Dirac/Weyl node, or refine.

---

### KI-12 — Refined-wedge star average is a different quadrature on non-orthogonal lattices

| | |
|---|---|
| Severity | low (documentation / expectation, not a wrong result) |
| Confidence | confirmed by reproducer |
| Location | `structure/es.py:_setCustomSymmetries` (docstring corrected in patch 13) |

The star average `M_sym(k) = (1/nsym) sum_S M(S k)` is **exact** on a regular
grid: the wedge sum reproduces the full-BZ sum to machine precision for any
lattice. The docstring further claimed it is exact on a *refined* wedge, via

```
sum_c (w_p/n^3) M_sym(k_c) = w_p * <M> over the union of the image cells
```

That collapse needs `{S k_c}` to be again the sub-cell centres of the image
cell `S.(parent cell)`, i.e. the `n x n x n` sub-tiling of a cell must be
invariant under the point group. True when the operations only permute and flip
the fractional axes (cubic / tetragonal / orthorhombic — KI-R3 was validated on
a cubic model). False for hexagonal, where a 60 degree rotation carries the
parallelogram patch of children onto a sheared one, so the union of the images
is not the refined reducible grid.

Graphene, 8x8x1 subdivided by 3, `sum_k w_k |v_xx|^2`:

```
refined reducible : 1.167906   (identical to the plain 24x24 grid)
refined wedge     : 1.169197   (0.11% apart)
converged 192x192 : 1.172925   (both ~0.4% away)
```

Both routes carry correct weights, exact group-averaged moments and conserved
total weight, and both converge to the same integral — they sample different
points. Nothing needs fixing in the code; what needed fixing was the claim.
The practical consequence is for *testing*: a refined wedge and the refined
reducible mesh it descends from must not be compared bitwise on a
non-orthogonal lattice. `testsuite/tests/test_refined_mesh_flags.py` asserts
machine precision on the coarse grid and only quadrature-level agreement after
refinement.

If bitwise agreement between the two routes is ever wanted, the subdivision
would have to be performed in a basis in which the point group acts by
permutation — or the children generated as the union over the star and then
re-reduced, which costs the whole point of working in the wedge.

---

### KI-13 — The refinement metric is necessary but not sufficient for transport convergence

| | |
|---|---|
| Severity | **high** (refinement can make the observable worse than its own coarse parent) |
| Confidence | confirmed by reproducer, on three lattices |
| Location | `structure/meshrefine.py:compute_error`, `select_refinement_panels` |

The sum-rule metric is a pure **energy-axis** criterion: it asks whether the set
of sampled band energies can integrate `df/da`. It never sees where in the
Brillouin zone those energies came from, so two meshes with identical energy
sets score identically however differently they sample k-space. Resolving the
kernel in energy is necessary for the transport integral to be accurate; it is
not sufficient.

This entry was previously rated *medium* and suspected of being an artefact of
a temperature/scale mismatch. Measurement (patch 20 session) shows the opposite:
the scale mismatch is real but secondary, and the blind spot is worse than
documented. Three independent symptoms, in increasing order of severity.

**1. The blind spot is specific to adaptive meshes.** On a *uniform* mesh the
energy axis and the k-space sampling are locked together, and the metric is a
good proxy — it tracks the transport error monotonically:

```
uniform cubic     nkp    metric(300K,10meV)   L11 error
8^3                35          0.2806            +145%
16^3              165          0.0326           +21.6%
24^3              455          0.0080            +8.2%
32^3              969          0.0050            +1.8%
44^3             2300          0.0043            +0.2%
56^3             4495          0.0043            +0.0%
```

Graphene behaves the same way (0.1866 / 0.0708 / 0.0245 / 0.0098 against
+142 / +22 / +3.3 / +0.6%). Adaptivity is what decouples the two, and only
then does the metric stop being informative. This is why
`scripts/kmesh_refinement.py:qualify_parent_mesh` trusts the metric for regular
grids and explicitly declines to trust it otherwise.

**2. At comparable metric, the transport error spans an order of magnitude.**
Graphene, `T = 300 K`, `gamma = 0.01 eV`, total `L11`, converged reference
`1.5432e+06`:

```
mesh                        nkp   metric(300K)      L11        error
uniform 48x48 irr           217      0.1866      3.7434e+06   +142%
adaptive <- 48 (iter 6)     441      0.0078      1.3525e+06   -12.3%
adaptive <- 48 (iter 7)     833      0.0037      1.3566e+06   -12.1%
adaptive <- 96             1201      0.0045      1.4294e+06    -7.4%
adaptive <- 192            3441      0.0044      1.5099e+06    -2.2%
adaptive <- 300            8011      0.0026      1.5281e+06    -1.0%
uniform 300x300 irr        7651      0.0098      1.5532e+06    +0.6%
```

The `adaptive <- 48` and `adaptive <- 192` meshes score 0.0078 and 0.0044 while
being 12.3% and 2.2% off. The uniform 300 mesh scores *worse* (0.0098) than
either and is 0.6% off. **The metric does not order the meshes by accuracy.**

Tightening the metric does not help: driving `error_tol` from 5e-3 to 5e-5 on
the graphene adaptive mesh saturates at 993 k-points and returns `1.3568e+06`,
i.e. it moves the observable by 0.02% while the deficit stays at 12%. The
residual is fixed by the *parent* mesh, not by refinement effort. The deviation
is likewise insensitive to `--defect_fraction` (0.9 / 0.99 / 0.999 give
-12.3 / -12.3 / -12.1%), so it is not a selection-tuning artefact.

**3. On an extended Fermi surface, refinement makes the observable WORSE than
its coarse parent.** Simple cubic at half filling, `T = 300 K`,
`gamma = 0.01 eV`, reference `2.0281e+07` (48^3 and 72^3 agree to 0.1%):

```
mesh                        nkp        L11         error
uniform 16^3 (the parent)   165     2.4655e+07    +21.6%
cascade   <- 16^3          1881     1.5493e+07    -23.6%
narrow-only <- 16^3        3675     1.4219e+07    -29.9%
uniform 44^3               2300     2.0332e+07     +0.2%
uniform 56^3               4495     2.0287e+07     +0.0%
```

Refining 165 -> 3675 k-points moves the result from +22% to -30%, and a plain
uniform mesh of comparable size is two orders of magnitude more accurate. The
same sign appears on the square lattice (16^2 parent +53.6%, narrow-only
-21.3%, cascade -16.4%).

The reason is structural: with a `(d-1)`-dimensional Fermi surface, band
energies near `mu` are dense *for free* — thousands of k-points already supply
them — so the energy-axis criterion is satisfied almost immediately while the
k-space sampling of the shell remains grossly inhomogeneous. The blind spot is
widest exactly where the Fermi surface is largest. Graphene's point node is the
easy case, and it is the case the method was tuned on.

**4. For an extended Fermi surface the narrow corner does not converge at
all.** Refining cubic to `(1 K, 1 meV)` from four different parents:

```
parent     nkp        L11(1K, 1meV)
8^3       7835         1.0966e+08
16^3      3675         5.0645e+07
24^3      4199         1.1800e+08
32^3      2477         7.7313e+07
```

A factor 2.3 of scatter, non-monotonic, no limit in sight. By contrast graphene
converges cleanly at its own corner (4.3182 / 4.3877 / 4.4279 / 4.4293 e+05
from parents 48 / 96 / 192 / 300, settling on 4.429e+05, with the 441-point
mesh only -2.5% off). So "the method works at the corner it was refined for"
holds for a point node and fails for a Fermi surface.

**What is not wrong.** The k-point weights and optical elements are exact on
every mesh tested: `sum_k w_k |v_xx|^2` is 2.000000 on all uniform meshes and
1.999995 on the adaptive ones. The error is entirely in how weight is
distributed *in energy*, not in the weights themselves. An earlier reading of
this entry blamed the weights; that reading is wrong.

**Mitigation shipped (patches 20, 21).** None of this is fixed, but it is no
longer silent:

- `--T_max` / `--gamma_max` run the refinement as a cascade over kernel widths,
  widest first, and the closing summary re-evaluates the final mesh at every
  rung so residual inadequacy at wide kernels is printed rather than implied.
  This improved graphene from -12.3% to -12.1% and cubic from -29.9% to -23.6%:
  visibility, not a cure.
- `qualify_parent_mesh` refuses a starting mesh that is too coarse at the widest
  kernel, with a concrete recommended `nk`, because that error is not
  recoverable by refining.
- `ltb inspect-mesh` / `lwann inspect-mesh` report density and axis proportions
  independently.

**Guidance until this is fixed.** Converge a *uniform* mesh at the widest kernel
of the planned sweep — cheap, and the metric is a valid proxy there — and use
refinement only to reach the narrow corner that no uniform mesh can afford (the
metric at 1 meV is still 0.267 on a 3001^2 graphene grid and falls as `1/n`).
Refinement is a low-temperature tool layered on a high-temperature-converged
parent, not a substitute for one. On an extended Fermi surface, prefer a finer
uniform mesh outright.

A weighted metric would be the natural remedy, but there is no exact sum rule to
anchor it: the weighted sum converges to `int DOS(e) K(e) de`, not to 1. The
options worth exploring are a self-consistency criterion (compare the transport
integral on the current mesh against the same integral on a once-more-refined
mesh) or Richardson extrapolation in the mesh density. Note that
`momentsDiagonal` is a single contiguous dataset, so a weighted criterion is
cheap to evaluate — it does not pay the per-k-point read penalty of KI-02.

---

### KI-14 — `lwann` has no check that the Hamiltonian respects the cell symmetry

| | |
|---|---|
| Severity | medium (wrong results, but only for hand-built low-symmetry input) |
| Confidence | confirmed by reading the code; the ltb analogue is confirmed by reproducer |
| Location | `structure/wannier.py`; cf. `structure/tb.py:_checkSymmetriesTightbinding` |

`ltb` verifies that the intra-orbital hoppings map onto themselves under the
point group spglib derives from the *structure*, and since patch 21 it raises
rather than warns when they do not (see KI-R17). `lwann` has no equivalent
check.

For genuine Wannier90 output this is a near-non-issue: the Hamiltonian was
derived from the same structure spglib is reading, so the two symmetries agree
by construction. The hazard is a hand-assembled or edited `_hr.dat`, or a model
written in Wannier format, whose hoppings are less symmetric than the cell. The
irreducible reduction would then fold k-points onto a star the band structure
does not possess. On the ltb side the same situation produced a 32% error in
`L11` on a test model before the check was escalated.

A hopping-pattern check like ltb's needs the orbital representation matrices
`D(g)`, which are real work for p/d manifolds. A cheaper representation-free
alternative is to sample a few hundred k-points and compare the *sorted*
eigenvalue spectrum at `k` and `g.k`: gauge-invariant, independent of orbital
count, and it tests exactly the property the reduction relies on,
`E_n(g k) = E_n(k)`. Cost scales as `norb^3`, so it would be reasonable as a
default up to a few orbitals and opt-in above. The same sampling could check
`M(g k) = g M(k) g^T` for the star average at negligible extra cost.

Partial mitigation exists: `structure/meshrefine.py:axis_anisotropy` warns when
symmetrising the measured axis anisotropy changes it by more than 10%, which is
the symptom of a stored point group that overstates the Hamiltonian's.

---

---

## Performance opportunities (identified, not implemented)

### PERF-02 — The generator re-reads and re-evaluates everything each iteration

| | |
|---|---|
| Expected gain | approaches the ratio (new points)/(total points) after a few iterations |
| Location | `structure/generators/*.py`, `structure/meshgenerator.py` |

Each refinement iteration constructs a fresh `TightBinding` /
`Wannier90Calculation`, re-parses the hopping file, and evaluates the model on
the **entire** mesh — including the k-points that were already evaluated in the
previous iteration and were not touched by the subdivision. After a few
iterations most of the mesh is unchanged, so most of the work is repeated.

**The mesh is purely additive.** With the enforced odd `refinement_factor`
the central child of a subdivided parent sits exactly on the parent, so every
parent coordinate survives; only its weight and `cell_delta` shrink. Verified
on a 4-parent test mesh: 4 parents in, 56 points out, **52 added, 0 removed**.
Since the per-k physics (energies, velocities, curvatures, moments) depends
only on the coordinate, every surviving row could be copied verbatim rather
than recomputed.

**What makes it non-trivial is bookkeeping, not physics** -- and *not*
symmetry: the star average at a k-point depends only on that k-point,
`symop` and `h(r)`, so a reused row stays valid on the symmetrised path too.
The obstacles are:

1. **Row indices are permuted.** `refine_kmesh` concatenates survivors and
   children and then calls `_merge_duplicates_with_tolerance`, which uses
   `np.unique` on rounded coordinates and therefore returns a
   lexicographically sorted mesh. A parent at row 0 came back at row 13 in the
   test above. An incremental generator needs an explicit
   `old_row -> new_row` map, which the refinement would have to return.
2. **The generator writes a fresh file from a fresh calculation object.**
   Reusing rows means reading the parent HDF5 back and splicing `energies`,
   `velocities`, `curvatures`, `moments` and `momentsBfield` into the new
   arrays -- and `momentsBfield` is stored as one dataset per k-point
   (KI-02), so the splice is `nkp` reads.
3. **Global quantities must still be recomputed** -- `_calcFermiLevel`, the
   gap detection, and the `opticalMoments[...,6:]` truncation test all look at
   the whole mesh. These are cheap, but they have to run after the splice, not
   before.
4. **`refinement_factor` must stay odd** for the additive property to hold.
   It is currently forced odd, but an incremental generator would depend on
   that invariant rather than merely benefit from it, so it should be asserted
   explicitly.

A fix therefore means extending the `MeshGenerator` protocol with an optional
`generate_incremental(new_points, new_weights, row_map, parent_hdf5,
output_path)` and having `refine_kmesh` return the row map. Worth doing now
that PERF-R1 is closed; the two touch the same code.

Cheaper intermediate step with most of the benefit for large Wannier models:
keep the calculation object alive across iterations (resetting the accumulator
lists) so the `_hr.dat` parse and the FT prefactor setup happen once.

### PERF-03 — A compiled / OpenMP kernel for the symmetrised path

| | |
|---|---|
| Expected gain | modest on top of PERF-R1 |
| Location | new |

`--openmp` currently only exports `OMP_NUM_THREADS` and friends before NumPy is
imported (`scripts/ltb/openmp_utils.py`), i.e. it is BLAS/LAPACK threading and
nothing else. There is no OpenMP anywhere in the refinement or in the
generators' own loops.

Now that PERF-R1 has vectorised the star expansion, the remaining hot spot is a
batch of small `eigh` calls plus a few `einsum` contractions, which BLAS
already threads. On the Fortran side, PERF-R2 removed the polygamma
bottleneck, so the same caution applies there: profile before assuming the
kernel is still where the time goes. A dedicated compiled kernel would only be worth it if
profiling after PERF-R1 still shows the per-k work dominating. Do not start
here.

### PERF-06 — `structure/tb.py` reducible branch is still unblocked

| | |
|---|---|
| Expected gain | memory only |
| Location | `structure/tb.py:_computeHk`, reducible branch |

The Wannier reducible branch was blocked in patch 05 and the symmetrising
branches of both files in patch 08, but the tight-binding *reducible* branch
still allocates `hk`, `hvk`, `hck` at full `(nkp, nbands, nbands[, 3/6])`.
For tight-binding models `nbands` is usually small enough that this does not
bite, which is why it was left alone; the blocking pattern from the other
three branches transfers directly if it ever does.

### PERF-07 — The chemical-potential search may dominate the polygamma budget

| | |
|---|---|
| Expected gain | unknown until measured; potentially large |
| Location | `src_linretrace/root.F90`, `occ_digamma_D/_Q` and the Ridders driver |

`occ_digamma` evaluates `psi_0` over `nband_max` (all bands, not the
`nbopt_min:nbopt_max` optical range) times `nkp`, **twice per Ridders
iteration**, with `maxiter` = 200 in the quad branch. `calc_polygamma`
evaluates `psi_1..psi_3` over `nbopt` times `nkp` **once**. The ratio of
evaluations is roughly

    2 * n_iter * nband_max / (3 * nbopt)

which for `nband_max = 20`, `nbopt = 10` and 20 iterations is ~27x more work
in the search than in the response.

**Measure before acting.** `niitact` is already printed per temperature step
next to `pot%MM`, and `timings(4)` already isolates polygamma evaluation from
`timings(5)`/`timings(6)`. Those two numbers settle whether this is the hot
spot at all.

If it is, the candidate is Newton rather than Ridders: the derivative is
analytic and already implemented as `polygamma2psi1` in `fermi.f90`,

    dn/dmu = (beta / 2 pi^2) * sum_k w_k * Z * Re psi_1(z)

so ~5 iterations instead of ~20, and `psi_range_*(z, 0, 1, psi)` returns both
orders from a single recurrence pass. Keep bisection as the fallback for the
`lT0 .and. muFermi` staircase case, where the existing comment correctly notes
that secant-type updates are ill-conditioned.

### PERF-04 — `openmp_utils` lives under `scripts/ltb/` but is generator-agnostic

| | |
|---|---|
| Expected gain | none (hygiene) |
| Location | `scripts/ltb/openmp_utils.py` |

`scripts/lwann/refine_kmesh.py` imports it from the `ltb` package. It should
move to `scripts/openmp_utils.py` with a re-export shim left behind, as part of
the planned refactor of the refinement into a fully input-agnostic module.

### PERF-05 — `detect_refinement_scale` probe grid is larger than it needs to be

**Superseded by patch 19 (KI-R14):** `detect_refinement_scale` no longer drives
the refinement; selection goes through `select_refinement_panels`, which needs
no probe grid at all. The function is retained only for callers that still
import it. Entry kept for provenance.

| | |
|---|---|
| Expected gain | small, already largely realised |
| Location | `structure/meshrefine.py` |

Patch 02 reduced the probe grid from 110k+10k to 12k+4k points, replaced the
element-wise `mpmath.polygamma` loop with a vectorised recurrence + asymptotic
series (~280x), and cached the analytic kernel across iterations. The function
now costs ~20–30 ms per iteration, so further tuning is not worth it unless the
probe grid is ever put back in the inner loop.

Retained here only to record the reasoning: the metric needs `max|defect|` and
its half-level set, both of which are resolved to well under a percent at the
reduced sizes.

---

## Resolved

Kept for provenance; do not renumber.

### KI-R1 — Boltzmann constant missing from the refinement metric

*Closed by patch 01.* `structure/meshrefine.py` used `beta = 1.0/temperature`
with the temperature documented and validated in Kelvin while all energies are
in eV. `beta` was therefore too small by `1/kB ~ 1.16e4`, i.e. the df/da kernel
was evaluated at a temperature four orders of magnitude too high. Fixed by
adding `kB_eV` to `structure/units.py` and routing both `digamma_im` and
`df_da` through a `_beta()` helper. See KI-05 for the remaining duplication.

### KI-R2 — Refined meshes reported `ndim = 0`

*Closed by patch 03/04.* Generators construct their calculation object with
`nkx = nky = nkz = 1`, so `_defineDimensions` reported `ndim = 0` and
`dims = [F,F,F]`. That propagated into `.unitcell/ndim`, where `linretrace`
warns about reduced dimensions and `postproc` takes its `ndim == 0` branch —
zeroing the off-diagonal Onsager elements and skipping the tensor inversion.
Every refined file produced before patch 03 carries wrong dimensionality
metadata and should be regenerated. See KI-07 for the residual caveat.

### KI-R3 — Irreducible meshes lost the optical-element symmetrisation

*Closed by patch 06.* `_computeHk` symmetrises the optical moments over the
star of each k-point only on its `irreducible` branch. The generators built
their calculation object with `irreducible=False`, so a refined wedge went
through the reducible branch: raw `v_i v_j` at wedge points, whose weighted sum
is not the BZ average. Cubic single-band model at 8x8x8:

```
irreducible run, symmetrised     : [0.2500  0.2500  0.2500]
same wedge via setCustomKmesh    : [0.2148  0.3203  0.2148]   <- wrong
reducible reference              : [0.2500  0.2500  0.2500]
```

Band energies are symmetry invariant, so the refinement error metric converged
normally and the failure was silent. Fixed by passing `.unitcell/symop` from
the coarse HDF5 into `setCustomKmesh`, which then selects the symmetrising
path. Validated: a refined *irreducible* mesh of 113 k-points and the
corresponding refined *reducible* mesh of 2280 k-points now agree to eight
decimal places.

**Follow-up:** the mechanism used to select that path — setting
`self.irreducible = True` in `_setCustomSymmetries` — also flipped
`.kmesh/irreducible` in the output file, which on the Fortran side means
something entirely different and broke the weights. See KI-R11.

### KI-R4 — Weight redistribution mixed inequivalent parents

*Closed by patch 01.* `refine_kmesh` pooled all removed parent weights and
redistributed them uniformly over all children of all parents. Total weight was
conserved, but weight was transferred between parents whenever their weights
differed — the generic case for an irreducible mesh, where the weight encodes
the star multiplicity. Replaced by a per-parent split `w_parent / n_children`.

### KI-R5 — Plain `lwann` crashed on every invocation

*Closed by patch 03.* The `--disorder` argparse block is commented out in the
public branch (the extension is not distributed) while `main()` still branches
on `args.disorder`, giving
`'Namespace' object has no attribute 'disorder'`. Fixed with
`parser.set_defaults(disorder=None)`, which is a no-op when the option is
re-enabled in a branch that ships the extension.

### KI-R6 — Log-spaced probe grid degenerated for `mu != 0`

*Closed by patch 02.* `detect_refinement_scale` built its mu-centred grid as
`geomspace(|mu|, |mu| + window)`. That is only genuinely logarithmic for
`mu ~ 0`: for `mu = 0.83, window = 0.1` the ratio is 1.12 and the "log" grid
degenerated into a second uniform grid, leaving the kernel peak unresolved.
Replaced by geometric spacing of the *offsets* from mu.

### KI-R7 — `structure.generators` was not a declared package

*Closed by patch 03.* `setup.py` did not list `structure.generators`, so
`LtbGenerator` was missing from any pip-installed build and `ltb refine` only
worked from a source checkout.

### KI-R8 — Wannier reducible branch stored a guaranteed-zero imaginary part

*Closed by patch 08.* For non-orthogonal lattices `opticalMoments` carries nine
components: six real `Re(v_a^* v_b)` and three imaginary parts of the
*off-diagonal* Cartesian pairs. `structure/tb.py` writes
`vk2[...,3:].imag` (components xy, xz, yz) in both of its branches, and so does
the irreducible branch of `structure/wannier.py` — but the reducible branch of
`structure/wannier.py` wrote `vel2[...,:3].imag`, i.e. the imaginary part of
`v_x^* v_x, v_y^* v_y, v_z^* v_z`, which is identically zero by construction.

Two consequences: the imaginary off-diagonal elements were silently lost, and
the truncation test that follows (`if not np.any(abs(opticalMoments[...,6:]) >
1e-6): truncate`) then *always* fired, dropping columns 6:9 entirely. Affected
`lwann` runs on non-orthogonal lattices only; orthogonal cells take the
three-component path and were never affected.

### KI-R9 — Wrong Bernoulli coefficient in the CERNLIB quad table

*Closed by patch 10/11.* `wpsipg_hp.F77` hard-codes all 80 series
coefficients `C(i,K)` as 30-digit literals. One is wrong:

```
C(13,0):  file  5.4827583333...Q+3
          exact 5.4827583333...Q+4      (B_26 / 26)      <- factor of ten
```

Every other entry checks out against `B_2i (2i+K-1)!/(2i)!`.

`K = 0` is `psi_0`, i.e. `occ_digamma`, i.e. the chemical-potential search.
The index 13 is only reached when `kmax >= 13`, so the **double** build
(`kmax_double = 7`) is unaffected and the **quad** build
(`kmax_quad = 16`) is. Measured worst relative error over K = 0..4 against
`mpmath` references, in real quad arithmetic:

```
wpsipghp (old)   3.41e-27
psi_fast (new)   9.81e-34        (quad eps = 1.93e-34)
```

with the discrepancy concentrated exactly where predicted -- `K = 0` at small
`|Im z|`:

```
z = (0.5, 0.0001)  K=0    old 2.8e-27    new 9.8e-34
z = (3.7, 0.02  )  K=0    old 3.4e-27    new 1.7e-34
```

So the quad branch has been running the chemical-potential search four digits
short of what it advertises.

`psi_fast.F90` does not transcribe the table. It stores 16 Bernoulli numbers
as exact integer numerator/denominator pairs and builds the K dependence from
exact small integers, because the factorial ratio is one
(`1/(2i)`, `1`, `2i+1`, `(2i+1)(2i+2)`, `(2i+1)(2i+2)(2i+3)`). Every number a
reader must check is now a Bernoulli numerator they can look up.

**Consequence for reproducibility:** quad results will not be bit-identical to
previous runs. The converged `mu` shifts at the ~1e-27 level -- unobservable,
but not zero.

### KI-R10 — Quad branch did not deliver quad precision

*Closed by patch 10/11.* Independently of KI-R9, the recurrence target and
the series order were mismatched. `wpsipg_hp` walks the argument up to
`Re V ~ 15` and then evaluates a 16-term Bernoulli series. Measured worst
case over K = 0..4 at `Re z = 1/2`, `Im z -> 0`:

```
walk to Re V ~ 15 (the old setting) : 4.3e-29
walk to |V|  = 22                   : 6.4e-37
walk to |V|  = 28                    : 2.1e-40
```

against a quad epsilon of 1.9e-34. The old code paid for 16 terms and then
stopped the walk too early to use them, losing ~5 digits precisely at
`E = mu` -- the one place where the quad branch exists to help.

`psi_fast` walks to `|V| = 28`, chosen as the smallest target at which an
argument that just fails the gate and one that just passes it both land within
4 x epsilon (`|V| = 28`, no walk, 16 terms: 2.2e-34). The double branch walks
to `|V| = 12` with a 9-term series, which is both more accurate and two steps
cheaper than the old `Re V ~ 15` with 7 terms.

### KI-R11 — `.kmesh/irreducible` conflated "wedge" with "regular grid", destroying refined-mesh weights

*Closed by patch 13 (Python) and patch 14 (Fortran).* **Severity: high — wrong
physical results, silently.**

`.kmesh/irreducible` has exactly one consumer in `linretrace`, and it is not
about symmetry at all. `input.f90` reads `/.kmesh/weights` only when the flag is
`.false.`, and `main.F90` otherwise reconstructs the weights in quad precision as

```fortran
nkred = kmesh%nkx*kmesh%nky*kmesh%nkz
kmesh%weightQ(ik) = kmesh%multiplicity(ik) * kmesh%weightsum / real(nkred,16)
```

The flag is therefore a promise that *the weights are reconstructible from the
grid dimensions*, which is true only for a regular Monkhorst-Pack grid.

The KI-R3 fix made `_setCustomSymmetries` set `self.irreducible = True` in
order to select the symmetrising branch of `_computeHk`. `h5output` writes that
same attribute, so every refined wedge went out flagged `irreducible = True`
while carrying the custom-mesh convention `nkx = nky = nkz = 1` and
`multiplicity = 1`. `linretrace` then discarded the (correct) `/.kmesh/weights`
and assigned `weightQ = weightsum` to **every** k-point.

Graphene, `templates/graphene.tbdata`, 12x12x1 coarse, three refinement steps,
`T = 300 K`, `gamma = 0.01 eV`, `L11` intra:

```
                          nkp        mu [eV]      L11_xx        L11_yy
reducible refined         192     2.5e-15      1.0939e+06    1.2040e+06
irreducible refined        43    -2.95        6.1484e+07    6.1484e+07   <- wrong
irreducible refined (fix)  43     2.5e-15     1.1490e+06    1.1490e+06
reducible 300x300       90000     6.3e-13     1.3772e+06    1.3772e+06
```

The total weight was `nkp * weightsum` instead of `weightsum`, so `L11` came out
a factor `nkp` too large and the chemical-potential search collapsed to the
bottom of the band. Reducible refined meshes were unaffected, which is why the
two routes disagreed.

**Python fix (patch 13).** The two meanings are now separate attributes:
`self.irreducible` keeps the grid-reconstructability meaning and is always
`False` for a custom mesh; the new `self.symmetrize` selects the star-averaging
code path. `h5output` writes `.kmesh/symmetrized` so that a continuation run
(`ltb refine refined_iter_N.hdf5 ...`) still picks up `symop`, and raises if
`irreducible` is set while `sum(multiplicity) != nkx*nky*nkz`.

**Fortran fix (patch 14).** `/.kmesh/weights` is now read unconditionally, and
the quad reconstruction is guarded by the same consistency condition
(`kmesh%uniformgrid`). Files written by the unpatched Python are therefore
rescued rather than silently mis-integrated. Verified: the legacy 43-point file
above reproduces the patched result bit for bit under the patched binary.

**Validation.** The symmetrised wedge equals the direction-average of the
reducible refinement of the same underlying mesh to 2e-11 relative:

```
refined reducible  (208 pts):  (L11_xx + L11_yy)/2 = 1181844.8306559924
refined irreducible (51 pts):   L11_xx = L11_yy    = 1181844.8306336943
```

Regular (ir)reducible grids are bit-identical before and after both patches.

### KI-R12 — `minimal_weight` was 1.0 for every custom mesh

*Closed by patch 14.* `input.f90` derived

```fortran
kmesh%minimal_weight = minval(kmesh%multiplicity) / real(kmesh%nkx*kmesh%nky*kmesh%nkz,8)
```

which for the custom-mesh convention collapses to `1/1 = 1.0`. It is the exit
criterion of the bisection in `find_mu_DFT` (`root.F90:128-129`), so with
`ElectronOccupation`, a `Bandgap` scissor, or scattering-file band shifts the
search returned after a single step, with an occupation deviation of up to one
full electron and a correspondingly arbitrary `mu_dft` and gap classification.
This affected **reducible** refined meshes too. Replaced by
`minval(inputweight) / weightsum`, which is algebraically identical for a
regular grid and correct for an adaptive one.

### KI-R13 — `.structure/weights` in the output file was uniform for adaptive meshes

*Closed by patch 14.* `output.F90` wrote `multiplicity / nkp`, which is uniform
for any custom mesh — exactly the adaptive weighting thrown away — and was
normalized inconsistently between reducible (sum = 1) and irreducible
(sum = nkx*nky*nkz/nkp) grids. `postproc/output.py` feeds this array to
`calcDOS`, so DOS/NOS plots made from the *output* file disagreed with the same
plots made from the *energy* file, which uses `.kmesh/weights`. Now written as
`kmesh%inputweight`, i.e. summing to `weightsum` in every case, so the NOS
saturates at the electron count. **Note:** this changes the normalization of
`.structure/weights` for regular grids as well; anything downstream that
re-normalized the old array should be checked.

### KI-R14 — Refinement metric was a signed sum and could rise on a finer mesh

*Closed by patch 19.* **Severity: medium — no wrong results, but the loop could
not converge and reported a false regression.**

`compute_error` returned `|1 - int df/da|` over the sampled band axis, i.e. the
absolute value of the SIGNED sum of per-panel trapezoid defects. The trapezoid
rule underestimates a concave integrand and overestimates a convex one, and
`df/da` is a peak: concave over its top, convex in both tails. The two
contributions therefore always carry opposite signs, in every system. When they
are of comparable size, an iterate whose peak defect cancels part of the tail
defect scores better than the strictly finer mesh that follows it.

Graphene 48x48x1 irreducible, `gamma = 1 meV`, `T = 1 K`:

```
nE        signed |sum d_i|    sum |d_i|
475            0.04480         0.04501
587            0.02301         0.04242   <- signed dips (cancellation)
1435           0.02966         0.03135   <- signed rebounds, L1 keeps falling
714051         0.00644         0.00665
```

The plateau detector read the 587 -> 1435 step as a **-28.89% regression** and
stopped. Refinement only adds k-points, so iteration 6 contained every point of
iteration 5 and was strictly the better sampling.

A second, independent defect made the target unreachable anyway.
`detect_refinement_scale` flagged energies by the POINTWISE interpolation defect
`|analytic - interpolated|` and converted their spread into a mu-centred radius.
The kernel amplitude is `1/(pi*gamma)` at mu against `gamma/(pi*eps^2)` in the
tails, a ratio of `(eps/gamma)^2`, so the peak always outranks the tails however
coarse the tails are. The window collapsed onto its `2*kernel_width` floor from
iteration 4 while the residual sat between 0.01 and 0.1 eV. Panel-resolved
budget at 714051 energies: peak `+0.00011`, tails `+0.00654`. The metric fell as
`N^-0.26`; with `--plateau_tol 0` the run reached only 0.0064 at 357217
k-points. Forcing the window open to its 0.1 eV ceiling did hit the target, but
needed 521905 k-points.

**Fix.** Both halves now use the same quantity, the per-panel quadrature defect
`d_i = trapezoid_i - exact_i` with `exact_i` from 5-point Gauss-Legendre:

* `compute_error` returns `sum|d_i| + |1 - sum exact_i|` (discretisation plus
  truncation of the sampled range). Monotone under refinement.
* `select_refinement_panels` marks panels by `|d_i|` until they account for
  `--defect_fraction` (0.9) of the total; `hotspot_mask_from_panels` selects the
  k-points bracketing them. Whole panels rather than a radius, so reaching a bad
  tail panel no longer drags in every k-point closer to mu.

Because `|sum d_i| <= sum|d_i|`, the metric is conservative: nothing previously
converged becomes unconverged. Runs that used to pass through cancellation now
need more iterations.

**Result.** Graphene, default settings, irreducible start:

```
start        iterations   final k-points   error
24x24x1           7            357         0.00112
48x48x1           6            425         0.00207
```

against *never converging* before. Controls unaffected: simple cubic 8^3
(gamma=1e-3, T=10K) converges in 3 iterations, cubic 16^3 (gamma=1e-2, T=300K)
in 1, square 16^2 in 2, and the `lwann refine` route in 7.

Rejected along the way: estimating `d_i` by the embedded rule
`(4/3)|trap_h - trap_h/2|`. It is 6x cheaper but collapses exactly where it is
needed, on a coarse mesh whose panels are wider than the kernel — 20.3 against a
true 39.7 on the graphene start. The kernel is closed-form, so Gauss-Legendre is
the right tool; 5 nodes is converged to four figures (0.006653 vs 0.006652 at 20
nodes) and costs 2.6 s at 714k panels.

### KI-R15 — Interband was documented as on by default and coded as off

*Closed by patch 16.* `config.f90:133` set `algo%lInterbandQuantities = .false.`
while `documentation/configspec:16` had always said `Interband = (default = T)`.
A hand-written config omitting the key silently produced an intraband-only run.
`lconfig` writes the key explicitly, so interactive users were unaffected.
Default is now `.true.`; the new `algo%lInterbandExplicit` (set from the
`bool_find` return) distinguishes an explicit request — missing full optical
elements remain a hard error — from the default being on, where main.F90 now
deactivates interband with a warning instead of aborting on `--intraonly` files.

### KI-R16 — `--energy_window` help text described a ceiling as the window

*Closed by patch 18.* `ltb refine --help` said "Energy window around mu for
hotspot detection. Default: 0.1", which reads as though 0.1 is the window in
force. It is a ceiling: the window used was
`min(energy_window, max(floor, reach))` and sat at its `2*kernel_width` floor
for most of a run. `lwann` already said "Ceiling for"; both now spell out the
formula and point at the per-iteration log line. (Patch 19 replaced the radius
mechanism entirely — see KI-R14 — but `energy_window` survives as the guard
against band-edge van Hove defects.)

### KI-R17 — Intra-orbital symmetry check had an inert orbital filter, and only warned

*Closed by patch 21.* Two defects in
`structure/tb.py:_checkSymmetriesTightbinding`.

First, the orbital filter never fired:

```python
band2_1  = int(self.tbdata[itb1][3]) - 1   # itb1 -- should be itb2
band2_2  = int(self.tbdata[itb1][4]) - 1   # itb1 -- should be itb2
if band1_1 != band2_1 or band2_2 != band2_2: continue   # name vs itself
```

Both identifiers were read from `itb1`, so `band2_1` was trivially equal to
`band1_1`, and `band2_2 != band2_2` is always false. A hopping belonging to
orbital pair `(a,a)` could therefore be matched against an entry from a
*different* orbital pair that happened to share an r-vector and a hopping
value, making the check too permissive on multi-orbital models — it could
report `True` on genuinely broken symmetry. Single-orbital models were
unaffected, which is why it went unnoticed.

Second, a failed check only warned and then built the irreducible mesh anyway.
On a test model whose hoppings broke the C4 symmetry of its cubic cell, that
produced `L11 = 2.886e+07` against the correct reducible-grid value of
`3.816e+07` — a 32% error that does not close with mesh density, because the
reduction folds each wedge point onto a star the band structure does not
possess. The check now raises on irreducible grids; `--red` remains the
documented escape and still works.

(For the record: the warning was working correctly and had been printing all
along. It was missed during the patch-20 investigation because the diagnostic
runs grepped stdout down to k-point counts.)

### PERF-R1 — Symmetrising branch was a Python loop over k-points

*Closed by patch 08.* See the patch notes: batching the star expansion plus
replacing the `einsum(optimize=True)` contractions with explicit BLAS calls
took a 512-k-point, 12-orbital, `nsym = 48` model from 36.5 s to 7.0 s.
Retained as a record of *why* the contractions are written out by hand: NumPy's
`einsum` re-plans per call, and its chosen paths made runtime and memory
non-monotonic in the batch size (33 s at `kblock=1`, 52 s at `kblock=2`, 7 s at
`kblock=8`, out-of-memory kill at `kblock=512`). Do not "simplify" them back
into `einsum`.

### PERF-R2 — Polygamma evaluation was several times slower than necessary

*Closed by patch 10/11.* Four independent problems in the CERNLIB routines, in
descending order of effect:

1. **The recurrence was entered on `|Re z|`, not `|z|`.** Every LinReTraCe
   argument has `Re z = 1/2 + beta*Gamma/2pi`, never near the threshold of 15,
   so the 14-step walk always ran -- but the step only increments the *real*
   part, so for a state with `|Im z| = 185` (a 1 eV band energy at 10 K) it
   moved `|V|` from 185.0007 to 185.6. Fourteen quad-precision complex
   divisions for nothing. Fraction of evaluations for which the walk was pure
   waste, `E - mu` over +-3 eV:

   | T | `|z| >= 28` |
   |---|---|
   | 1 K | 99.5 % |
   | 10 K | 94.6 % |
   | 100 K | 45.9 % |
   | 300 K | ~0 % |

2. **Three separate calls for K = 1,2,3 on the same argument.**
   `calc_polygamma` looped `do ipg = 1,3`, and each call redid the reduction
   and recomputed `1/V**(K+1)` from scratch: 45 complex divisions where 15
   divisions plus a few multiplies suffice. In quad this is disproportionate
   -- measured on this hardware, a complex division is ~625 ns against ~184 ns
   for a multiply, both software-emulated by libquadmath.

3. **The expansion order was a compile-time constant** (`-D kmax=`), which
   also controlled *which* `DATA` statements existed through a nest of `#if`
   blocks. `psi_fast` always defines all 16 coefficients and picks the Horner
   start index at run time from `|V|`; at T = 1 K most states need 5-6 terms
   rather than 16. No preprocessor symbols remain.

4. **`occ_digamma` uses only `aimag(psi_0)`**, and `psi_0` is the one order
   needing a complex logarithm -- ~2046 ns in quad, against ~663 ns for the
   `atan2` that yields the imaginary part alone. `psi0_imag_*` replaces it.

Measured in real Fortran quad, `E - mu` over +-3 eV, `Gamma` = 1 meV:

| | psi_1..psi_3 old -> new | psi_0 old -> new |
|---|---|---|
| 1 K | 41.3 -> 7.4 us (5.5x) | 13.2 -> 3.6 us (3.7x) |
| 10 K | 43.3 -> 10.2 us (4.3x) | 14.6 -> 4.8 us (3.0x) |
| 100 K | 44.3 -> 22.7 us (2.0x) | 14.4 -> 9.4 us (1.5x) |
| 300 K | 44.0 -> 37.1 us (1.2x) | 14.6 -> 14.3 us (1.0x) |

The gain grows as temperature falls and as the relevant `E - mu` window
widens -- which is where the quad branch matters, and the direction the
topological transport kernels are heading, since terms proportional to the
Fermi function have no `df/dw` cutoff.

