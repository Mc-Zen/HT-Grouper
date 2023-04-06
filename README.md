# Hardware-tailored grouping for Pauli measurements




## Dependencies

- Gurobi

## Setup
- Check that you have Gurobi installed and registered a valid license
- Download or clone this repository
- Run CMake at repository level
  - During configuration, you may choose the path for the Gurobi binaries by setting the variable `GUROBI_PATH` (defaults to `C:/gurobi1000/win64`)
  - Generate and open project
- The grouper algorithm can be found and run from the `ht-grouper` sub-project
