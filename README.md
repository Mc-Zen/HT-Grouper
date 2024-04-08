# HT-Grouper 

This repository is maintained by Kyano Levi and Daniel Miller [(AG Eisert, FU Berlin)](https://www.physik.fu-berlin.de/en/einrichtungen/ag/ag-eisert/people/index.html).

## Grouping Pauli operators into commuting sets that admit hardware-tailored readout circuits

In the preprint [Hardware-Tailored Diagonalization Circuits](https://doi.org/10.48550/arXiv.2203.03646), Daniel Miller _et al._ developed a theoretical framework for the construction of hardware-tailored readout circuits.
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

- Gurobi >= 10.0.0
- Qiskit >= 0.36.0
- Numpy >= 1.22

 The code has only been tested on Windows. 


## Setup

This is a quick setup guide for users that are somewhat familar with C++ development. Check out the [full installation guide](docs/installation-guide.md) for detailed step-by-step information on how to setup this code. 

A recent C++ compiler version is needed, supporting the C++20 standard. We highly recommend using Visual Studio 2019 or 2022 (not to be confused with Visual Studio Code).

- Check that you have Gurobi installed and registered a valid license.
- Download or clone this repository.
- Run CMake at repository level, configure, generate and open project. 
- The grouper algorithm can be found and run from the `grouper` sub-project.


## Workflow Example

We have prepared a jupyter notebook to guide you through the workflow: [data/WorkflowExample.ipynb](data/WorkflowExample.ipynb)


## Performance

Below, we showcase a benchmark of HT-Grouper carried out on a standard laptop (Windows Intel i7, 8th generation, 8 threads).
The code was executed in release mode (fast option).

The benchmarked Hamiltonians describe $n$-atomic hydrogen chains in STO-3G basis mapped to $2n$ qubits via the Bravyi-Kitaev fermi-to-qubit mapper. 
The readout circuits were tailored to linear connectivity.
We applied HT-Grouper with three hyperparameter choices: 

- (bright green) all subgraphs --> exponential runtime 
- (medium bright green) 10k subgraphs --> polynomial runtime
- (medium dark green) 1k subgraphs --> polynomial runtime
- (dark green) 100 subgraphs --> polynomial runtime

As expected, the [estimated shot reduction](https://doi.org/10.22331/q-2021-01-20-385) lies between the grouping result of [Sorted Insertion](https://doi.org/10.22331/q-2021-01-20-385) with general commutativity (red), which requires all-to-all connectivity, and Sorted Insertion with qubit-wise commutativity (blue), which only requires single-qubit Clifford gates at the readout stage.

Here, we did not fully explore the tradeoff between runtime and quality:
For example, for 24 qubits and 7151 Paulis (10 hydrogen atoms), the HT-grouper with 100 random subgraphs required less than 35 minutes to terminate. However, the resulting estimated shot reduction is 12.4 (still better than blue), which can likely be improved by re-running the HT-Grouper with an increased number of subgraphs. 

![grafik](https://github.com/Mc-Zen/HT-Grouper/assets/129524538/ab8ce32a-1227-40c5-94d4-7bbe5ab2d1b9)


