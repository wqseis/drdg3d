﻿# DRDG3D [![DOI](https://zenodo.org/badge/593392070.svg)](https://zenodo.org/badge/latestdoi/593392070)

For 3D dynamic rupture (DR) modelling, using a nodal discontinuous Galerkin (DG) finite element method [[1]](#1).

## Installation

Here we show how to install the software in MacOS. For the installation on Linux system, please refer to [manual_drdg3d.pdf](https://github.com/wqseis/drdg3d/blob/main/doc/manual_drdg3d.pdf).

Before the compilation of the software, make sure you have [make](https://www.gnu.org/software/make/), [gcc](https://gcc.gnu.org/), [openmpi](https://www.open-mpi.org/), [netcdf-fortran](https://docs.unidata.ucar.edu/netcdf-fortran/current/), [gmsh](https://gmsh.info/), [meshio](https://pypi.org/project/meshio/) and [metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) installed.

 ```bash
> brew install make
> brew install gcc
> brew install open-mpi
> brew install netcdf
> brew install netcdf-fortran
> pip install gmsh
> pip install meshio
> pip install metis
```

### Compilation

After installing requirements, the program is installed as follows.

```bash
> git clone https://github.com/wqseis/drdg3d.git
> cd drdg3d
> mkdir -p bin obj
> make
```

Three executable files will be generated in the bin  directory: **exe_solver**, **exe_get_neigh** and **exe_part_mesh**.

Before running the examples, set the path first:
```bash
> bash scripts/setenv.sh
```

## Usage

We take tpv5 as an example:

```bash
> cd examples/tpv5
```
1. Generate the mesh:

```bash
> bash genmesh.sh
```

2. Run the preprocess in MATLAB:

    **run_preprocess.m**

3. Run the program
```bash
> bash go.sh
```

### Plot the results

Use **draw_fault_snap.m** to draw the snapshots of fault slip rate (or slip, stresses, etc.). There are other scripts such as **draw_fault_seismo.m**, **draw_grdsurf_snap.m**, **draw_grdsurf_seismo.m**, please refer to [manual_drdg3d.pdf](https://github.com/wqseis/drdg3d/blob/main/doc/manual_drdg3d.pdf) for more details.

## Reference

<a id="1">[1]</a> Wenqiang Zhang, Yajing Liu, Xiaofei Chen. A Mixed-Flux-Based Nodal Discontinuous Galerkin Method for 3D Dynamic Rupture Modeling. *_ESS Open Archive_.* November 01, 2022. DOI: [10.1002/essoar.10512657.1](https://www.doi.org/10.1002/essoar.10512657.1)

<a id="2">[2]</a> Zhang, W., Liu, Y., & Chen, X. (2023). A mixed-flux-based nodal discontinuous Galerkin method for 3D dynamic rupture modeling. Journal of Geophysical Research: Solid Earth, 128, e2022JB025817. [https://doi.org/10.1029/2022JB025817](https://doi.org/10.1029/2022JB025817)

