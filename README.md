# LinReTraCe
<!-- [<img src="https://github.com/LinReTraCe/LinReTraCe/blob/release/documentation/logo.png?raw=true" width="196" height="196">](https://sites.google.com/view/tomczak-group/projects/linretrace) -->

The [Linear Response Transport Centre](https://sites.google.com/view/tomczak-group/projects/linretrace) (LinReTraCe) is a package for the simulation of transport properties driven by carriers with finite lifetimes. The underlying theory, described in [PRB:105.085139](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.105.085139) ([arxiv:2112.07604](https://arxiv.org/abs/2112.07604)), establishes a comprehensive and thermodynamically consistent phenomenology capable of reproducing qualitatively correct full temperature profiles in metals as well as semiconductors. A comprehensive code documentation including implementation details, benchmarks, and test cases is described in [SciPostPhysCodeb.16](https://scipost.org/10.21468/SciPostPhysCodeb.16)
and a step-by-step installation guide is available in the [guide](https://github.com/LinReTraCe/LinReTraCe/blob/release/documentation/guide.pdf).

The code package provides several interfaces to common electronic structure codes, including [Wien2K](http://susi.theochem.tuwien.ac.at), [VASP](https://vasp.at), as well as maximally localized Wannier functions from [Wannier90](http://www.wannier.org). The DFT input can also be supplied with the band interpolation scheme of [BoltzTraP2](https://gitlab.com/sousaw/BoltzTraP2).

Moreover, we provide an interface to create general tight-binding models, as well as a Python3 interface where source-agnostic data can be supplied via specifically shaped arrays.

## Prerequisites

At its core, LinReTraCe is a highly efficient and scalable MPI parallalized Fortran code for the calculation of accurate, artefact-free transport coefficients of, both, realistic electronic structures as well as models and with high precision down to lowest temperatures. All the surrounding interfaces and tools are written in modern Python3.
In order to obtain all required and optional packages for the running of the pre- and postprocessing at its full functionality, simply install the dependencies via [pip](https://pypi.org/project/pip/). To avoid package incompatibilities with global installations, we highly recommend the usage of local environments. In the root folder of linretrace use

`python -m venv pylrtc`

`source pylrtc/bin/activate`

`pip install -e .`

If you want to include the BoltzTraP2 features also use

`pip install cmake packaging pyfftw boltztrap2`

To return to your standard environment, simply use

`deactivate`

## Compilation

The Fortran part of LinReTraCe requires a full [HDF5](https://www.hdfgroup.org/solutions/hdf5/) installation (`version >=1.12.1`) whose underlying HDF5 library calls are handled with an [HDF5 wrapper](https://github.com/linretrace/hdf5_wrapper) written by one of us. At the wrapper's page, an installation guide and test code for the HDF5 library is provided, if needed (Please note that one needs to register at the HDF5 site to download the tar ball). To maximize the scalabilty of the code, we recommend making use of the MPI implementation.

To obtain the LinReTraCe source, clone this repository:

`git clone https://github.com/linretrace/linretrace`

To compile the code, a special `make_include` file needs to be saved in the main linretrace folder. An examplary configuration that enables MPI looks as follows
```
FC       = mpiifort       # Main Fortran compiler
FCDG     = ifort          # Fortran Compiler for the Special functions
FFLAGS   = -O3
FPPFLAGS = -DMPI
HDF5     = -I/opt/hdf5-1.13.1_icc/include
HDF5    += -L/opt/hdf5-1.13.1_icc/lib -lhdf5_fortran -lhdf5hl_fortran
```
A single core configuration looks as
```
FC       = gfortran
FCDG     = gfortran
FFLAGS   = -O3
HDF5     = -I/opt/hdf5-1.13.1_gcc/include
HDF5    += -L/opt/hdf5-1.13.1_gcc/lib -lhdf5_fortran -lhdf5hl_fortran
```
The configuration can be checked via `make validate` and the compilation is done with `make linretrace`, creating the `bin` subfolder in which the binary `linretrace` will be moved into.
An overview of validation, compilation, and testsuite is displayed via `make`.

## Workflow
![LinReTraCe workflow](https://github.com/LinReTraCe/LinReTraCe/blob/release/documentation/flowchart.png?raw=true "LinReTraCe workflow")

### Energy file
The center point of LinReTraCe's input is the energy file, where all the necessary band energies as well as optical elements (among other auxiliary data) are stored. In order to prepare this file one of the following interfaces can be used:

**WIEN2k** and **VASP** are interfaced with `ldft`:

`ldft <wien2k folder> --optic --output wien2k_structure.hdf5`

`ldft <vasp folder> --interp --output vasp_structure.hdf5`

In the WIEN2k example we make use of the dipole matrix elements from the optic subpackage (`x optic`) whereas in the VASP example we interface the BoltzTraP2 band interpolation scheme via `--interp`. There the Peierls approximation is used based on optical elements consisting of band velocities and curvatures. For a full description, see the code publication and `ldft --help`.

**Wannier90** is interfaced with `lwann`:

`lwann <wannier90 folder> --output wannier90_structure.hdf5`

Please note that we provide the possibility to expand the reducible grid with `--kmesh nkx nky nkz` as well as sub-interface generic WIEN2k momemtum grids if the required `case.struct` and `case.klist` files are provided:
`lwann <wannier90 folder> --wien2k --output wannier90_wien2k.hdf5`

General **tight-binding models** are created with `ltb`. Simply provide a text file with the corresponding hopping parameters, atomic positions and lattice vectors and execute:

`ltb tb_file nkx nky nkz charge --output model.hdf5`

The tight binding file has a Wannier90 inspired format:
```
begin hopping
#  a1 a2 a3    orb1 orb2  hopping.real [hopping.imag]
   0  0  0     1    1     0.3  # on site energy
  +1  0  0     1    1     1.0  # nearest neighbor hopping
   0 +1  0     1    1     1.0
  -1  0  0     1    1     1.0
   0 -1  0     1    1     1.0
end hopping

begin atoms
#  sort rx ry rz
   1    0  0  0     # fractional lattice vector coordinates
end atoms

begin real_lattice
#  x   y   z
   5   0   0       # a1 lattice vector in units of Angstroem
   0   5   0       # a2
   0   0   1       # a3
end real_lattice
```
Please note that we use the convention employed by the strongly correlated electron systems community, where a positive hopping leads a reduction of the energy, i.e.

$$ H(\mathbf{k}) = -\sum_\mathbf{R} e^{i\mathbf{k}\cdot\mathbf{R}} (1-2\delta_{\mathbf{R},\mathbf{0}}\delta_{l,l'}) H_{ll'}(\mathbf{R}) $$

If instead one wants to provide energies and optical elements that cannot be created with the above tools we provide a **generic interface** `linterface` that contains the class `StructureFromArrays` that supports the load-in of the necessary data (multiplicity, energies, optical elements, ...).

### Adaptive k-mesh refinement

Transport coefficients are integrals over a kernel of half-width
$\max(k_BT,\Gamma)$ centred on the chemical potential. When that width is
small compared with the band dispersion, a uniform mesh dense enough to
resolve it is wasteful almost everywhere else. `ltb refine` and `lwann refine`
grow a non-uniform mesh instead, subdividing only where the integration
actually fails:

`ltb refine coarse.hdf5 tb_file <mu> <gamma_min> <T_min> [options]`

`lwann refine coarse.hdf5 <wannier90 folder> <mu> <gamma_min> <T_min> [options]`

`gamma_min` and `T_min` are the *narrow corner* of the sweep you intend to run
afterwards: the smallest scattering rate and the lowest temperature, which
together set the narrowest kernel the mesh will have to resolve. Note that
$\Gamma$ floors the width — at $T = 10$ K, $k_BT = 0.86$ meV, so a run at
$\Gamma = 10$ meV is Lorentzian-dominated at both ends of a 10–300 K sweep and
the kernel width spans only a factor 2.6.

Each iteration diagonalises the model on the current mesh, measures how well
the sampled band energies integrate $\partial f/\partial \varepsilon$, marks
the energy intervals carrying the quadrature error, and subdivides the
k-points that bracket them by `--refinement_factor` (3) per active axis. Cell
widths are carried in `/.kmesh/cell_deltas` so that weights stay exact and a
run can be continued from any `refined_iter_N.hdf5`.

Irreducible coarse meshes are supported and are the cheaper route: the wedge
is kept, and the optical elements are averaged over the star of every k-point.
Refined meshes are written with `.kmesh/irreducible = False` (they are not
regular grids, so their weights come from `.kmesh/weights`) plus a
`.kmesh/symmetrized` marker so continuation runs keep symmetrising.

#### Read this before relying on a refined mesh

**Refinement cannot repair a starting mesh that is too coarse at the widest
kernel of your sweep.** The accuracy there is set by the parent and is not
recoverable afterwards, because at the narrow corner those energy intervals
carry no selectable quadrature defect. Measured on graphene, tightening
`--error_tol` by a factor 100 moved the 300 K result by 0.02%, while changing
the parent from 48x48 to 300x300 moved it from -12.3% to -1.0%.

**On an extended Fermi surface, refinement can make the observable worse than
its own coarse parent.** Simple cubic at half filling, 300 K, $\Gamma = 10$ meV,
against a converged reference of `2.0281e+07`: the 16^3 parent (165 k-points)
gives +21.6%, refining it to 3675 k-points gives **-29.9%**, and a plain
uniform 44^3 mesh (2300 k-points) gives +0.2%. The metric is an energy-axis
criterion, and with a two-dimensional Fermi surface the energy axis near $\mu$
is dense for free while the k-space sampling of the shell stays inhomogeneous.

So the recommended workflow is:

1. Converge a **uniform** mesh at the *widest* kernel of the planned sweep.
   This is cheap and the metric is a valid proxy on a regular grid, where it
   tracks the transport error monotonically.
2. Refine from that mesh only to reach the **narrow** corner, which is the one
   no uniform mesh can afford (the metric at 1 meV is still 0.267 on a
   3001x3001 graphene grid and falls only as $1/n$).

Refinement is a low-temperature tool layered on a high-temperature-converged
parent, not a substitute for one. Where the Fermi surface is extended, a finer
uniform mesh is often simply better. See KI-13 in `KNOWN_ISSUES.md` for the
full measurements.

#### Declaring the wide corner: the cascade

`--T_max` and `--gamma_max` declare the *widest* corner of the sweep. The
refinement then runs once per rung of a logarithmic ladder of kernel widths,
**widest first**, so that the wide stage places k-points across the whole
$\pm W_\mathrm{max}$ shell while it is still cheap and those points become the
parents the narrow stages subdivide. Taking the maximum over the ladder inside
a single loop would not do this: on a coarse mesh the narrow-kernel defect
dominates, so the maximum is always the narrow entry.

```
ltb refine coarse.hdf5 model.tb 0.0 1e-3 10 --T_max 300 --gamma_max 0.01
```

Without these flags the behaviour is unchanged and the mesh is certified at the
narrow corner only. With them, the starting mesh is first **qualified** at the
wide corner and the run refuses with a recommended `nk` if it is too coarse
(`--parent_tol`, default `--error_tol`; `--skip_parent_check` to bypass). The
closing summary re-evaluates the final mesh at every rung, so residual
inadequacy at wide kernels is printed rather than implied.

Other options: `--ladder_ratio` (3, matching `--refinement_factor`),
`--ladder_max_stages` (6), `--stage_tol_exponent` (0, i.e. the same tolerance
at every width — this grades the mesh self-similarly, which is what the
metric's $(h/W)^2$ scaling asks for), and `--max_output_size` (e.g. `4GB`),
which predicts the next file from the current bytes-per-k-point *before* the
generator runs and stops, keeping the last mesh that fitted.

#### Inspecting a mesh

`ltb inspect-mesh` and `lwann inspect-mesh` report what a mesh can resolve at
one kernel width, answering two questions that need two different criteria:

```
ltb inspect-mesh mesh.hdf5 --T 300 --gamma 0.01 --target 5e-3
```

*Density* is judged by the sum-rule metric, and a recommended `nk` is given
when it falls short. *Axis proportions* are judged separately, by the
band-diagonal velocity tensor, because the metric is blind to them — it ranks a
32x32 mesh above a 64x8 one that is 250x more accurate at the same cost.

The proportions matter more than they look. On an orthorhombic model with
$t_x/t_y = 5$, at matched cost: 1:1 gives +12.3%, 2:1 +2.1%, 4.6:1 +0.50%, and
8:1 +0.004%. Note this is **not** the usual DFT rule $n_i \propto |b_i|$: in a
tight-binding model the hoppings are given per lattice vector, so stretching a
cell axis does not flatten the band along it. What matters is the energy
variation per fractional step, which is what `inspect-mesh` measures. The
reported ratio is a lower bound and is rounded up in the suggestion; the
accuracy plateau is broad, so over-shooting is close to free.

Exit code 2 means "below target", for scripting.

See `documentation/adaptive_kmesh.tex` for the algorithm and `ltb refine
--help` for all options.

### Scattering File
Arbitrary (momentum, band, and spin dependent) scattering rates (quasi particle weights and energy shifts) are supported through a custom HDF5 scattering file. The workflow to create a custom file follows
- Copy `lscat_template` from the installation folder into your working direction.
- Insert the linretrace folder into the system path.
- Reference to correct energy file.
- Define calculation axis ($\mu$ or $T$-scan).
- Define scattering rates as numpy array.
Optionally: Define quasi particle weights and/or band shifts.
- Execute script to generate LRTC scattering file.

Simplistic temperature dependencies (polynomial) on the other hand can be generated via the config file.
### Config File
LinReTraCe itself is configured via a free format configuration file. `lconfig` provides a minimal starting point through interactive questioning. Its main purpose is to define which quantities should be calculated at which precision in addition to mode-specific sub configurations. The temperature mode, e.g. further supports impurity levels, impurity bands, and homogeneous doping. For a full description see `documentation/configspec`.
### LinReTraCe
Running your `LinReTraCe` installation is done either via
`mpirun -np <cores> bin/linretrace config.lrtc`
or
`bin/linretrace config.lrtc`
where `config.lrtc` is the config file from the previous subsection.

### Output File
The generated output is an HDF5 file, whose tree structure is documented in [paper]. To interface this file in an effortless way we provide `lprint`, capable of plotting/printing all the available transport coefficients.

`lprint <LRTCoutput file> list`
lists all the available datasets, which then can be retrieved by providing the corresponding key and, optionally, directional arguments:

`lprint -p <LRTCoutput file> s-intra xx yy`
plots, e.g., the xx and yy components of the intra-band Seebeck tensor

`lprint -p <LRTCoutput file> rh-intra xyz`
plots the Hall coefficient in the xy-plane (magnetic-field in z-direction)
## License
This project is licensed under the GNU General Public License v3 which can be found in **LICENSE**.


## Authors
Matthias Pickem\*, Emanuele Maggio, Jan M. Tomczak\*

\* Corresponding authors:

matthias [dot] pickem [at] gmail [dot] com

tomczak [dot] jm [at] gmail [dot] com

## Acknowledgements
[LinReTraCe](https://sites.google.com/view/tomczak-group/projects/linretrace) was funded by the Austrian Science Fund (FWF) through project [P 30213](https://pf.fwf.ac.at/de/wissenschaft-konkret/project-finder/40827).
