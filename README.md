# HT-Grouper 

This repository is maintained by Kyano Levi and Daniel Miller [(AG Eisert, FU Berlin)](https://www.physik.fu-berlin.de/en/einrichtungen/ag/ag-eisert/people/index.html).

## Grouping Pauli operators into commuting sets that admit hardware-tailored readout circuits

In the preprint [Hardware-Tailored Diagonalization Circuits](https://doi.org/10.48550/arXiv.2203.03646) by Daniel Miller _et al._, 
a theoretical framework for the construction of hardware-tailored readout circuits was developed.
The original implementation of the accompanying algorithms is not publically available.
However, the working principle of the grouping algorithms are precisely described in [Sec. II](https://doi.org/10.48550/arXiv.2203.03646)  of the supplementary material.
This enabled us to reimplement an open-source version of one of the algorithms. 
Specifically, in this GitHub repository, we implement a modified version of the the algorithm based on the "numerical solver" described in [Sec. II D ](https://doi.org/10.48550/arXiv.2203.03646)  of the supplementary material.

This repository contains a C++ implementation of the algorithm that relies on the Mixed Integer Quadratically Constrained Program (MIQCP).
Our implementation leverages [Gurobi](https://www.gurobi.com/downloads/gurobi-software/), 
a commercially available software with [special offers for academics](https://www.gurobi.com/academia/academic-program-and-licenses/).


Furthermore, this project includes C++ code for dealing with matrices, undirected graphs, Pauli operators, Clifford circuits and more. 

## License

This library is distributed under the [MIT License][license].

If you want to support work like this, please cite our paper:
[arXiv.2203.03646](https://doi.org/10.48550/arXiv.2203.03646)

[license]: https://github.com/Mc-Zen/HT-Grouper/blob/master/LICENSE.txt


## Dependencies

- Gurobi 10.0.0 (or more recent)
- Qiskit 

 The code has only been tested on Windows. 


## Setup

This is a quick setup guide for users that are somewhat familar with C++ development. Check out the [full installation guide](docs/installation-guide.md) for detailed step-by-step information on how to setup this code. 

A recent C++ compiler version is needed, supporting the C++20 standard. We highly recommend using Visual Studio 2019 or 2022 (not to be confused with Visual Studio Code).

- Check that you have Gurobi installed and registered a valid license.
- Download or clone this repository.
- Run CMake at repository level, configure, generate and open project. 
- The grouper algorithm can be found and run from the `grouper` sub-project.


## Code example 1

[TODO]

show code that finds HT circuits for {P1, P2}, {P1,P3}, {P2,P3} but fails for {P1, P2, P3}:
![grafik](https://github.com/Mc-Zen/HT-Grouper/assets/129524538/aa53f136-46a7-499f-8a2b-480df6df35f4)


## Code example 2

[TODO]

show code that reconstructs the Pauli grouping in Tab IX for the Hamiltonian in Tab VII.

Also apply sorted-insertion-qwc and print R^_HT/R^_TPB.
