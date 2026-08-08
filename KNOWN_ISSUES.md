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

Last updated: 2026-08-08 (patch 09).

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
already threads. A dedicated compiled kernel would only be worth it if
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

### PERF-04 — `openmp_utils` lives under `scripts/ltb/` but is generator-agnostic

| | |
|---|---|
| Expected gain | none (hygiene) |
| Location | `scripts/ltb/openmp_utils.py` |

`scripts/lwann/refine_kmesh.py` imports it from the `ltb` package. It should
move to `scripts/openmp_utils.py` with a re-export shim left behind, as part of
the planned refactor of the refinement into a fully input-agnostic module.

### PERF-05 — `detect_refinement_scale` probe grid is larger than it needs to be

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

### PERF-R1 — Symmetrising branch was a Python loop over k-points

*Closed by patch 08.* See the patch notes: batching the star expansion plus
replacing the `einsum(optimize=True)` contractions with explicit BLAS calls
took a 512-k-point, 12-orbital, `nsym = 48` model from 36.5 s to 7.0 s.
Retained as a record of *why* the contractions are written out by hand: NumPy's
`einsum` re-plans per call, and its chosen paths made runtime and memory
non-monotonic in the batch size (33 s at `kblock=1`, 52 s at `kblock=2`, 7 s at
`kblock=8`, out-of-memory kill at `kblock=512`). Do not "simplify" them back
into `einsum`.

