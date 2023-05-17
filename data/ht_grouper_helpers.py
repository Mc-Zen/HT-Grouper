import json
from typing import Dict, List
from qiskit import QuantumCircuit
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import InverseCancellation
from qiskit.circuit.library import HGate


h_gate_canceller = PassManager([InverseCancellation([HGate()])])


def read_hamiltonian(filename: str) -> Dict[str, float]:
    """
    Read a hamiltonian from a JSON file containing Pauli strings
    as keys and real-valued coefficients as values. 

    Pauli strings are read as 
       `"XYZ"` -> X on qubit 0, Y on qubit 1, Z on qubit 2

    Returns a dictionary of the same format. 
    """
    with open(filename) as file:
        return json.load(file)


def write_hamiltonian(filename: str, hamiltonian: Dict[str, float], invert_pauli_order: bool = False):
    """
    Write a hamiltonian in form of a dictionary with Pauli strings
    as keys and real-valued coefficients as values to a JSON file. 

    Pauli strings should read as 
       `"XYZ"` -> X on qubit 0, Y on qubit 1, Z on qubit 2
    If they are read in the inverse order, set `invert_pauli_order` to `True`
    """
    if invert_pauli_order:
        hamiltonian = {key[::-1]: value for key, value in hamiltonian.items()}
    with open(filename, "w") as file:
        json.dump(hamiltonian, file, indent=2)


def read_grouping(filename: str) -> List[dict]:
    """
    Read a Pauli grouping from a JSON file of the following format:
    ```
    {
      "grouping": [
        {
          "operators": ["III","IIZ","ZZZ"],
          "edges": [],
          "cliffords": ["H","H","H"]
        },
        {
          "operators": ["XYZ","YZI","IIX"],
          "edges": [[0,1],[1,2]],
          "cliffords": ["S", "HS", "I"]
        }
      ]
    }
    ```

    A grouping is a list of groups, each containing a list of Pauli 
    operators, a list of edges (stored as 2-element lists) representing
    the graphs (or CZ gates) and a list of Clifford gates that need to be 
    applied to each qubit in the readout circuit. 

    Valid Clifford gates are "I", "H", "S", "SH", "HS" and "HSH". 

    Pauli strings are read as 
       `"XYZ"` -> X on qubit 0, Y on qubit 1, Z on qubit 2


    Parameters
    ----------
    filename : str
        _description_

    Returns
    -------
    List[dict]
        _description_
    """
    with open(filename) as file:
        result = json.load(file)
        return result["grouping"]


def generate_readout_circuits(grouping: List[dict]) -> List[QuantumCircuit]:
    """
    Generate readout circuits from a Pauli grouping specified in the format
    as described in the documentation of :func:``read_grouping()``
    """
    circuits = []
    clifford_map = {"H": QuantumCircuit.h, "S": QuantumCircuit.s, "I": lambda *params: None}
    for group in grouping:
        num_qubits = len(group["operators"][0])
        edges = group["edges"]
        cliffords = group["cliffords"]
        circuit = QuantumCircuit(num_qubits)
        for qubit, clifford in enumerate(cliffords):
            for gate in clifford[::-1]:
                clifford_map[gate](circuit, qubit)

        for q1, q2 in edges:
            circuit.cz(q1, q2)
        circuit.h(range(num_qubits))

        circuits.append(h_gate_canceller.run(circuit))
    return circuits


def R_hat(grouping: List[dict], hamiltonian: Dict[str, float]) -> float:
    """
    Compute estimated shot reduction compared to single Pauli measurements. 

                       ∑_i^N ∑_j^{m_i} |a_ij|    
         \hat{R} = (----------------------------)²  
                         ∑_i^N √(∑_j^{m_i} |a_ij|²) 

        as defined in https://doi.org/10.22331/q-2021-01-20-385


    Parameters
    ----------
    grouping : List[dict]
        Pauli grouping, a list of dictionaries each containing a list of 
        Pauli strings for the key "operators", i.e. grouping[0]["groupings"]
        needs to return a list of Pauli strings. 

    hamiltonian : Dict[str, float]
        A hamiltonian specification as a dictionary with Pauli strings as keys 
        and coefficients as values

    Returns
    -------
    float
        Estimated shot reduction R_hat
    """
    numerator = 0
    denominator = 0
    identity = "I" * len(next(iter(hamiltonian)))

    for group in grouping:
        term = 0
        for pauli in group["operators"]:
            if pauli == identity:
                continue
            absolute = abs(hamiltonian[pauli])
            numerator += absolute
            term += absolute * absolute
        denominator += term ** .5
    return (numerator / denominator) ** 2


hamiltonian = read_hamiltonian("hamiltonians/H_chain.json")

grouping = read_grouping("groupings/ham2_grouping.json")
grouping = read_grouping("groupings/H_chain_grouping.json")

# circuits = generate_readout_circuits(grouping)
# for circuit in circuits:
#     print(circuit)

print(R_hat(grouping, hamiltonian))
