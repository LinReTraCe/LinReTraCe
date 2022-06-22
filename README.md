# LinReTraCe
<!-- [<img src="https://github.com/LinReTraCe/LinReTraCe/blob/release/documentation/logo.png?raw=true" width="196" height="196">](https://sites.google.com/view/tomczak-group/projects/linretrace) -->

The [Linear Response Transport Centre](https://sites.google.com/view/tomczak-group/projects/linretrace) (LinReTraCe) is a package for the simulation of transport properties driven by carriers with finite lifetimes. The underlying theory, described in [PRB:105.085139](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.105.085139) ([arxiv:2112.07604](https://arxiv.org/abs/2112.07604)), establishes a comprehensive and thermodynamically consistent phenomenology capable of reproducing qualitatively correct full temperature profiles in metals as well as semiconductors. A comprehensive code documentation including implementation details, benchmarks, and test cases is described in [arxiv:2206.06097](https://arxiv.org/abs/2206.06097)
and a step-by-step installation guide is available in the [userguide](https://github.com/LinReTraCe/LinReTraCe/blob/release/documentation/userguide.pdf).

The code package provides several interfaces to common electronic structure codes, including [Wien2K](http://susi.theochem.tuwien.ac.at), [VASP](https://vasp.at), as well as maximally localized Wannier functions from [Wannier90](http://www.wannier.org). The DFT input can also be supplied with the band interpolation scheme of [BoltzTraP2](https://gitlab.com/sousaw/BoltzTraP2).

Moreover, we provide an interface to create general tight-binding models, as well as a Python3 interface where source-agnostic data can be supplied via specifically shaped arrays.

## Prerequisites

At its core, LinReTraCe is a highly efficient and scalable MPI parallalized Fortran code for the calculation of accurate, artefact-free transport coefficients of, both, realistic electronic structures as well as models and with high precision down to lowest temperatures. All the surrounding interfaces and tools are written in modern Python3.
In order to obtain all required and optional packages for the running of the pre- and postprocessing at its full functionality, simply install the dependencies with either [pip](https://pypi.org/project/pip/)

`pip install matplotlib h5py numpy scipy ase spglib`

`pip install boltztrap2`

or use one of many other popular Python package managers like [anaconda3](https://www.anaconda.com).

The Python interface of LinReTraCe can then either bei executed directly from the source folder or installed globally via

`python setup.py install`

## Compilation

The Fortran part of LinReTraCe requires a full [HDF5](https://www.hdfgroup.org/solutions/hdf5/) installation (`version >=1.12.1`) whose underlying HDF5 library calls are handled with an [HDF5 wrapper](https://github.com/linretrace/hdf5_wrapper) written by one of us. At the wrapper's page, an installation guide and test code for the HDF5 library is provided, if needed. To maximize the scalabilty of the code, we recommend making use of the MPI implementation.

To obtain the LinReTraCe source, clone this repository:

`git clone https://github.com/linretrace/linretrace.git`

To compile the code, a special `make_config` file needs to be saved in the main linretrace folder. An examplary configuration that enables MPI looks as follows
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
The compilation is done with `make`, creating the `bin` subfolder in which the binary `linretrace` will be moved into.

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
