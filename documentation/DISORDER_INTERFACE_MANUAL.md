# Ab-initio disorder interface for LinReTraCe — workflow manual

`lwann --disorder` turns a two-point disorder archive (`deltaH.h5`, produced
by the disorder-Wannier pipeline) into LinReTraCe input: momentum- and
band-resolved on-shell scattering rates, virtual-crystal bandshifts, and
transition matrix elements dressed with the resummed static disorder
(diffuson) ladder. The LinReTraCe Fortran core is not modified in any way.

Physics and derivations: `documentation/disorder_interface.pdf`.
This manual covers operation: requirements, workflow, every flag, every
output dataset, the diagnostics, and troubleshooting.

---

## 1. Requirements

- The LinReTraCe repository (this tree) with `lwann` and
  `structure/wannier_disorder.py`.
- The **disorder-Wannier pipeline modules** (`dwannier.py`,
  `disorder_selfenergy.py`, `disorder_transport.py`). Point to their
  directory with the environment variable `DWANNIER_PATH` or the flag
  `--dwannier-path`.
- A disorder archive `deltaH.h5` (pipeline format >= 1.1). It stores the
  two-point disorder Hamiltonian **per configuration** (real-space blocks
  plus weights); ensemble means and fluctuations are formed on load by
  `DisorderArchive` — see the notes, Sec. 2, for why per-configuration
  storage is the right (lossless) choice and how the interface exports the
  averaged correlator for analysis.
- A Wannier90 folder for the **primitive** cell (`case.nnkp`, `case.win`,
  `case.wout`, `case_hr.dat`), the same Wannier basis the archive was built
  from. The interface verifies H(k) consistency at runtime and reports the
  deviation.

## 2. The workflow in three commands

```bash
export DWANNIER_PATH=/path/to/dwannier_modules

# (1) plain energy file on a q_p-commensurate, Gamma-centred k-mesh
lwann w90folder --kmesh 6 6 6 -o lrtc-wann.hdf5

# (2) disorder input: scattering file + dressed energy file + lambda file
lwann w90folder --kmesh 6 6 6 -o lrtc-wann.hdf5 \
    --disorder deltaH.h5 \
    --disorder-engine scba \
    --disorder-temps 100 500 9 --disorder-temps-log \
    --vertex-eval onshell \
    --scat-output lrtc-dis-scat.hdf5 \
    --dressed-energy-output lrtc-wann-dressed.hdf5 \
    --lambda-output lrtc-lambda.hdf5

# (3) run the UNTOUCHED core on the pair (dressed energies, disorder rates):
#     in the linretrace config point the energy file to
#     lrtc-wann-dressed.hdf5 and the scattering file to lrtc-dis-scat.hdf5
linretrace config.lrtc
```

If step (2) is invoked without an existing energy file, the plain run of
step (1) is performed automatically first.

**k-mesh rule.** Every archive transfer must land on the mesh:
`q_p * (nkx, nky, nkz)` must be integer for all `p` (half-integer transfer
components need even divisions), and the mesh must be Gamma-centred
(unshifted). Violations abort with an explanatory message.

**Combined with DMFT** (`--dmft w2d_T1.hdf5 w2d_T2.hdf5 ...` together with
`--disorder`): the DMFT temperature axis takes precedence over
`--disorder-temps`; static parts add before the band rotation; the rails of
rates and ladder carry the total width Z·(Gamma_dmft + Gamma_dis); the rung
stays disorder-only; lambda becomes temperature dependent (computed per T).
The dressed energy file carries the first temperature's vertex (a warning
reminds you; per-T dressed files are future work).

## 3. All flags

| flag | default | meaning |
|---|---|---|
| `--disorder HDF5` | — | the archive; enables everything below |
| `--disorder-engine {born,scba,ata}` | `scba` | on-shell rate engine: single Born pass / self-consistent Born (the conserving partner of the ladder) / average T-matrix (per-configuration multiple scattering, bound-state trends) |
| `--disorder-temps TMIN TMAX NT` | `100 500 5` | temperature axis in K (elastic rates are T-independent and replicated) |
| `--disorder-temps-log` (alias `--disorder-tlog`) | off | logarithmic T mesh, as in the linretrace config; ignored when `--dmft` is present |
| `--no-disorder-vca` | off | do **not** fold the archive's virtual-crystal Sigma^(1)(k) into rotation and bandshifts |
| `--no-disorder-vertex` | off | rates/shifts only; skip the ladder and the dressed energy file |
| `--rate-eval {onshell,fermi}` | `onshell` | evaluation energy of the on-shell rates: at eps=a_kn (quasiparticle width; activated regime) or eps=0 (tail weight at mu; metals, low-T saturation). Keep matched with `--vertex-eval` (warning otherwise) |
| `--vertex-channels {ladder,cooperon,both}` | `ladder` | which vertex channels dress the matrix elements: diffuson ladder, maximally crossed (Cooperon / weak localisation), or both. Cooperon requires an inelastic dephasing rate (`--gamma-add` or `--dmft`) and triggers Q-mesh / supercell / ensemble adequacy warnings |
| `--gamma-add C0 C2` | `0 0` | phenomenological INELASTIC rate Gamma_inel(T)=C0+C2*T^2 [eV, eV/K^2], added consistently to rate rails, ladder rails, the written scatrate (per T step) and acting as the Cooperon dephasing cutoff. Makes the dressing T-dependent |
| `--dressed-per-t` | off | with a T-dependent dressing (combined mode or `--gamma-add`): one dressed energy file per temperature (suffix `_T####`); run linretrace once per file |
| `--vertex-eval {onshell,fermi}` | `onshell` | evaluation energy of the static ladder: `onshell` = at eps = a_kn per state (activated regime of semiconductors; retains particle-hole asymmetry for L12/Seebeck); `fermi` = at eps = 0, one solve per orbit (metals; low-T tail-dominated saturation regime). **Run both**: disagreement flags the bimodal crossover where a single static dressing is strained (see notes Sec. 6) |
| `--coherence-window X` | `6.0` | dress inter-band elements for \|a_n − a_m\| < X·(Gt_n + Gt_m); outside, bare elements are kept |
| `--disorder-svd-cut X` | `1e-2` | relative truncated-SVD cut deflating the near-unit conserved-density (diffusion) direction of the DC ladder |
| `--scat-output F` | `lrtc-disorder-scat.hdf5` | scattering HDF5 name |
| `--dressed-energy-output F` | `lrtc-wann-dressed.hdf5` | dressed energy HDF5 name |
| `--lambda-output F` | off | write the lambda analysis HDF5 |
| `--dwannier-path DIR` | `$DWANNIER_PATH` | location of the pipeline modules |

Shared with the rest of lwann: `--kmesh`, `-o/--output`, `--energy-hdf5`
(reuse an existing energy file with a non-default name), `--nocorrection`
(skip the inter-orbital Peierls correction), `--dmft ...` and all DMFT fit
flags.

## 4. Outputs — where to find which quantity

### 4.1 Scattering HDF5 (`--scat-output`)
FullScattering layout (identical to the DMFT interface); per spin prefix
(`/`, or `/up/`,`/dn/`):

```
.quantities/tempAxis, betaAxis, muAxis, nT, ...
step/000001/scatrate    (nkp, nbands)   bare rate [eV]  (core multiplies by qpweight on read)
step/000001/qpweight    (nkp, nbands)   Z (=1 disorder-only)
step/000001/bandshift   (nkp, nbands)   from Sigma^(1)(k) [+ DMFT shifts in combined mode]
... one step per temperature
```

### 4.2 Dressed energy HDF5 (`--dressed-energy-output`)
Byte-copy of the plain energy file with `momentsDiagonal` and the
in-coherence-window rows of `kPoint/*/moments` replaced by the one-sided
ladder dressing `Re[conj(v^a) Lambda^b]`. Attributes `disorder_dressed`
(eval mode) and `disorder_archive`. **This is the file to feed to
linretrace** together with the scattering file.

### 4.3 Lambda analysis HDF5 (`--lambda-output`)

| dataset | shape | content |
|---|---|---|
| `lambdaDiagonal` | (nkp, nbands, ndir) | dressed/bare ratios; guarded to 1 where the bare moment is negligible (band extrema, v→0 — the dressed values remain exact). Leading (nT, …) axis in combined mode |
| `vertexDiagonal` | (nkp, nbands, 3) complex | Lambda^alpha_nn itself |
| `kPoint/NNNNNNNNNN/lambdaPairs` | (npairs, 8) | rows (n, m, Re Lambda^{x,y,z}, Im Lambda^{x,y,z}) of dressed inter-band pairs |
| `correlatorW2` | (nkp, N, nb, nb) | ensemble-averaged squared correlator sum_c w_c \|dV_band\|^2 in the band basis on the actual mesh; [k, p, n, m] couples (k, n) to (k − q_p, m). Exported for analysis (the archive itself stores per-configuration data — lossless for all higher moments) |
| `bandEnergies`, `railWidths` | (nkp, nb) | a_kn and total rail widths used in the solves |
| attrs | — | archive path, `vertex_eval`, temperature axis |

## 5. Runtime diagnostics and how to react

Every run prints (grep-able):

- **`archive/lwann H0(k) consistency`** — max deviation between the archive
  reference Hamiltonian and the lwann Fourier transform on sampled k. Should
  be ~1e-12 eV or better; larger values mean the Wannier runs behind the two
  inputs differ (energies/moments then come from lwann, disorder quantities
  from the archive — reconcile the inputs).
- **`correlator parity asymmetry`** — violation of inversion symmetry of the
  ensemble-averaged correlator. Above 5% the ensemble is not self-averaging:
  the parity protection of the current channel is weakened; treat lambda
  qualitatively and add configurations.
- **`Gamma_dis in [...]`** and the **level-spacing warning** — if the median
  rail width is below the median level spacing of the mesh, on-shell
  quantities are k-mesh limited: densify `--kmesh`.
- **`vertex[...]: lambda min/max/mean, dropped modes, cond, source leak`** —
  the deflation report. Leak > 1e-3 triggers a warning: the current source
  genuinely drives the (dropped) diffusion mode — broken parity protection;
  lambda near the affected energies is unreliable (rates and shifts are not
  affected). Remedy: more configurations.
- **`lambda ratio guard`** — count of (k, band, dir) entries where the bare
  moment is negligible and lambda was set to 1 for reporting; the dressed
  matrix elements themselves are exact everywhere.

## 6. Validation

```bash
DWANNIER_PATH=/path/to/dwannier \
DWANNIER_ARCHIVE=/path/to/deltaH.h5 \
DISORDER_W90DIR=w90folder \
python3 testsuite/tests/test_disorder_interface.py
```

Checks: orbit ladder solver == the machine-validated pipeline solver
(`disorder_transport.solve_vertex`) at ~1e-14 on the ab-initio correlator;
point-like correlator gives lambda = 1 exactly (parity); forward-scattering
correlator matches the analytic tau_tr/tau to 1e-10; end-to-end run in both
`--vertex-eval` modes with dressed == lambda·bare at machine precision.

## 7. Inelastic rates and the Cooperon channel

Three routes for scattering of non-disorder origin (e-e, e-ph — inelastic):

1. **`--gamma-add C0 C2`** (preferred with vertex corrections): enters the
   rate rails, the ladder rails, the written scatrate (per T step), and
   provides the Cooperon dephasing cutoff — fully consistent.
2. **Config `ScatteringCoefficients` / `QuasiParticleCoefficients`
   together with `ScatteringFile`** — enabled by the small core extension
   in this patch series (config.f90 + input.f90): the polynomial
   Gamma_add(T) = sum_i c_i T^(i-1) (T in Kelvin) is ADDED per temperature
   step on top of the file's scatrate, and the config Z is combined with
   the file's qpweight via 1/Z_tot = 1/Z_file + 1/Z_cfg − 1, before the
   Gamma·Z convention. Kernels consistent; the dressed elements were built
   without this rate (fine for Gamma_inel << Gamma_dis).
3. **`ScatteringOffset`** (unmodified core): constant, added then ×Z.
   Same consistency caveat as (2).

**Cooperon workflow** (`--vertex-channels cooperon|both`): supply an
inelastic rate first (the run refuses otherwise — an unregularised Q→0
pole; elastic disorder does not dephase). The run prints the length
scales (elastic mean free path l_el, dephasing length L_phi, supercell
size, Q-mesh resolution) and warns when: dQ_max > 1/L_phi (mesh, not
physics, cuts the WL logarithm — increase `--kmesh` to N_i ≥ |b_i|·L_phi);
l_el exceeds the supercell size (coherent disorder replicas — increase
det S so the supercell linear size ≳ l_el); the ensemble is small or
parity-asymmetric (pole region amplifies ensemble noise). The crossed
contribution is T-dependent through tau_phi(T): use `--dressed-per-t`
and run linretrace per temperature. The Q loop parallelises with the
pipeline's tile_map (`DWANN_NUM_WORKERS`, `DWANN_BACKEND=mpi`) — same
rationale as the full scheme.

## 8. Core extension shipped in this series

`src_linretrace/config.f90` + `input.f90`: with a `ScatteringFile`, the
config keys `ScatteringCoefficients` / `QuasiParticleCoefficients` are now
parsed and applied ADDITIVELY per temperature step (see route 2 above).
The change mirrors existing idioms and touches nothing else; note the
export used here did not contain the `.F90` preprocessor sources, so
please run one compile + regression on your side (the text-mode
polynomial convention in main.F90 should match T in Kelvin,
coefficients c_i in eV/K^(i-1)).

## 9. Troubleshooting

| symptom | cause / remedy |
|---|---|
| `k-mesh ... not commensurate with the archive transfers` | choose `--kmesh` so that every q_p·mesh is integer (even meshes for half-integer transfers) |
| `requires an unshifted k-mesh` | remove the Monkhorst-Pack shift |
| `archive n_p != Wannier90 nproj` | the archive and the w90 folder use different Wannier bases |
| `energy file ... inconsistent with this lwann run` | regenerate the energy file on the same `--kmesh` (or pass it via `--energy-hdf5`) |
| ImportError about dwannier modules | set `DWANNIER_PATH` / `--dwannier-path` |
| huge lambda spread + leak warning | non-self-averaging ensemble (few configurations): rates/shifts fine, vertex qualitative — enlarge the ensemble |
| `Cooperon ... unregularised` error | supply `--gamma-add C0 C2` or `--dmft` (inelastic dephasing) |
| `Q-mesh too coarse for the dephasing length` | increase `--kmesh` to the printed N_i ≥ \|b_i\|·L_phi |
| `mean free path exceeds the supercell` | larger supercell (det S) so L_sc ≳ l_el |
| lambda differs strongly between `onshell` and `fermi` | bimodal-weight crossover regime; benchmark against the full `disorder_transport` machinery |
| overwrite prompt in scripts | the FullScattering writer asks interactively; from Python use `computeDisorder(..., interactive=False)` |

## 10. Limitations / future work

- Disorder-current end terms (need the ∂ₖΔH data): available in the full
  pipeline (`disorder_transport`), percent-level at DC on the test systems.
- Per-temperature dressed energy files in combined DMFT mode.
- B-field (Hall) matrix elements are copied through undressed (the
  three-current vertex has its own structure; use the full machinery for
  vertex-corrected Hall).
