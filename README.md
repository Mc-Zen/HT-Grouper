# Hardware-tailored grouping for Pauli measurements

This repository contains an implementation of the algorithm presented in [this paper](https://doi.org/10.48550/arXiv.2203.03646) for prioritized grouping of Pauli operators into commuting subsets that can each be measured simultaneously by a readout circuit that is tailored to a certain quantum hardware connectivity. 

Also, this project features a variety of C++ libraries for dealing with matrices, undirected graphs, Pauli operators, clifford circuits and more. 


## Dependencies

- Gurobi

## Setup
- Check that you have Gurobi installed and registered a valid license
- Download or clone this repository
- Run CMake at repository level
  - During configuration, you may choose the path for the Gurobi binaries by setting the variable `GUROBI_PATH` (defaults to `C:/gurobi1000/win64`)
  - Generate and open project
- The grouper algorithm can be found and run from the `ht-grouper` sub-project
