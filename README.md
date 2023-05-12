# Hardware-Tailoring Pauli Grouper

This repository is maintained by Kyano Levi and Daniel Miller [(AG Eisert, FU Berlin)](https://www.physik.fu-berlin.de/en/einrichtungen/ag/ag-eisert/people/index.html).

## Grouping Pauli operators into commuting sets that admit hardware-tailored readout circuits

In the preprint [Hardware-Tailored Diagonalization Circuits](https://doi.org/10.48550/arXiv.2203.03646) by Daniel Miller et al., 
a theoretical framework for the construction of hardware-tailored readout circuits was developed.
The original code of the presented algorithms is not publically available, however, the algorithms are precisely described in Sec. II of the supplementary material.
Here, we present an open-source version of the algorithm that is referred to as ``numerical solver'' in [arXiv.2203.03646](https://doi.org/10.48550/arXiv.2203.03646).

This repository contains a C++ implementation the algorithm that relies on the Mixed Integer Quadradically Constrained Program (MIQCP).
Our implementation leverages [Gurobi](https://www.gurobi.com/downloads/gurobi-software/), 
a commercially available software with [special offers for academics](https://www.gurobi.com/academia/academic-program-and-licenses/).


Furthemore, this project includes C++ code for dealing with matrices, undirected graphs, Pauli operators, Clifford circuits and more. 


## Dependencies

- Gurobi


## Setup

A recent C++ compiler version is needed, supporting the C++20 standard. Using Visual Studio 2019 or 2022 is recommended. 

- Check that you have Gurobi installed and registered a valid license
- Download or clone this repository
- Run CMake at repository level
  - During configuration, you may choose the path for the Gurobi binaries by setting the variable `GUROBI_PATH` (defaults to `C:/gurobi1000/win64`)
  - Generate and open project
- The grouper algorithm can be found and run from the `ht-grouper` sub-project


## Code example

[TODO]
