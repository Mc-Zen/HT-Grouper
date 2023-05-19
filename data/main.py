from ht_grouper_helpers import *


def main():





    hamiltonian = read_hamiltonian("hamiltonians/H_chain.json")

    grouping = read_grouping("groupings/ham2_grouping.json")
    grouping = read_grouping("groupings/H_chain_grouping.json")

    # circuits = generate_readout_circuits(grouping)
    # for circuit in circuits:
    #     print(circuit)

    print(R_hat(grouping, hamiltonian))

    qc = QuantumCircuit(8)
    qc.h(range(8))

    exp = HamiltonianExperiment(qc, grouping, hamiltonian)
    result = exp.run(20000)
    expectation_values = exp.evaluate(result)
    print(expectation_values)
    
main()