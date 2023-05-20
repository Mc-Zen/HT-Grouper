import numpy as np
import qiskit
from qiskit.providers import Job
from qiskit.quantum_info import Pauli
from qiskit.result import Result
import json
from typing import Dict, List, Optional, Sequence, Union
from qiskit import QuantumCircuit
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import InverseCancellation
from qiskit.circuit.library import HGate


h_gate_canceller = PassManager([InverseCancellation([HGate()])])


def read_hamiltonian_from_json(filename: str) -> Dict[str, float]:
    """
    Read a hamiltonian from a JSON file containing Pauli strings
    as keys and real-valued coefficients as values. 

    Pauli strings are read as 
       `"XYZ"` -> X on qubit 0, Y on qubit 1, Z on qubit 2

    Returns a dictionary of the same format. 
    """
    with open(filename) as file:
        return json.load(file)


def write_hamiltonian_to_json(filename: str, hamiltonian: Dict[str, float], invert_pauli_order: bool = False):
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


def read_grouping_from_json(filename: str) -> List[dict]:
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


Bitstring = np.int64


class BinaryResult:
    """
    Class for storing a single quantum circuit result with result
    signature and count. 

    The result signature is encoded in a binary string in little-endian
    (least-significant bit represents the zeroth quantum bit result). 

    """

    __slots__ = ("bitstring", "count")

    def __init__(self, bitstring: Union[Bitstring, int], count: int):
        self.bitstring = Bitstring(bitstring)
        self.count = count

    def __eq__(self, other) -> bool:
        return self.bitstring == other.bitstring and self.count == other.count

    def __str__(self, num_qubits=0) -> str:
        return f"{self.bitstring:0>{num_qubits}b}: {self.count}"

    def __repr__(self) -> str:
        return f"BinaryResult({self.bitstring}, {self.count})"


class CircuitResult:
    """
    Class encapsuling results from one quantum circuit, similar to 
    :class:`qiskit.result.ExperimentResultData` but the keys are not stored
    using hexadecimal strings but binary bitstrings which allows for
    some performance improvements when evaluating readout results. . 

    """

    __slots__ = ("results", "num_qubits")

    def __init__(self, counts: Dict[str, int], qubits: Optional[Sequence[int]] = None):
        """Create a CircuitResult from a dictionary as returned by `qiskit.result.get_counts()`.
        The keys are little-endian bitstrings, i.e. the zeroth register is at the rightmost 
        position in the string.

        E.g., ``{"001": 5, "010": 4}`` to denote a result where the outcome with only the zeroth 
        qubit being 1 occured 5 times and the outcome of the first qubit begin 1 occured 5 times. 

        Optionally, the results on only a selection of qubits can be extracted (marginalization).  
        Example: 

        >>> CircuitResult({"101": 9, "001": 4}, [0, 2])

        will result in storing ``{"11": 9, "01": 4}``. 

        Parameters
        ----------
        counts : Dict[str, int]
            Outcome/count dictionary
        qubits : Optional[Sequence[int]], optional
            If specified, only the given qubits are extracted. 
        """
        self.results: List[BinaryResult] = []

        if qubits is None:
            if len(counts) != 0:
                self.num_qubits = len(next(iter(counts)).replace(" ", ""))
            for key, value in counts.items():
                self.results.append(BinaryResult(Bitstring(int(key.replace(" ", ""), 2)), value))

        else:
            self.num_qubits = len(qubits)
            for key, value in counts.items():
                key = key.replace(" ", "")  # might contain spaces to separate registers
                key = "".join(key[index] for index in qubits)
                self.results.append(BinaryResult(Bitstring(int(key, 2)), value))

    def __str__(self) -> str:
        return "CircuitResult(" + ", ".join(result.__str__(self.num_qubits) for result in self.results) + ")"


class HamiltonianExperiment:

    def __init__(self, preparation_circuit: QuantumCircuit, grouping: List[dict], hamiltonian: Dict[str, float]):
        self.preparation_circuit = preparation_circuit
        self.grouping = grouping
        self.hamiltonian = hamiltonian
        self.num_qubits = len(grouping[0]["operators"][0])
        assert self.num_qubits == preparation_circuit.num_qubits, "Number of qubits do not match for preparation circuit and hamiltonian"

    def get_readout_circuits(self) -> List[QuantumCircuit]:
        if not hasattr(self, "readout_circuits"):
            self.readout_circuits = generate_readout_circuits(self.grouping)
        return self.readout_circuits

    def get_circuits(self) -> List[QuantumCircuit]:
        if not hasattr(self, "circuits"):
            readout_circuits = self.get_readout_circuits()
            self.circuits = []
            for readout_circuit in readout_circuits:
                qc = self.preparation_circuit.compose(readout_circuit)
                qc.measure_all()
                self.circuits.append(qc)
        return self.circuits

    def simulate(self, shots):
        circuits = self.get_circuits()
        backend_sim = qiskit.Aer.get_backend("qasm_simulator")
        job = qiskit.execute(circuits, backend=backend_sim, shots=shots)
        return job

    def evaluate(self, job: Job) -> Dict[Pauli, float]:
        all_counts = job.result().get_counts()
        expectation_values: Dict[Pauli, float] = {}

        for group, readout_circuit, counts in zip(self.grouping, self.get_readout_circuits(), all_counts):
            circuit_result = CircuitResult(counts)
            for pauli_string in group["operators"]:
                pauli = Pauli(pauli_string[::-1])
                pauli_z = pauli.evolve(readout_circuit, frame="s")
                assert (pauli_z.phase == 2 or pauli_z.phase == 0) and not pauli_z.x.any()

                expectation_value = _compute_expectation_value(circuit_result, create_bitstring_from_nparray(pauli_z.z))
                if pauli_z.phase == 2:
                    expectation_value *= -1
                expectation_values[pauli] = expectation_value

        expectation_values[Pauli("I" * self.num_qubits)] = 1.
        return expectation_values


def create_bitstring_from_nparray(arr: np.ndarray) -> Bitstring:
    bitstring = Bitstring(0)
    for i in range(len(arr)):
        if arr[i]:
            bitstring |= (1 << i)
    return bitstring


def _compute_expectation_value(circuit_result: CircuitResult, s: Bitstring):
    expectation_value: int = 0
    total_count: int = 0
    for result in circuit_result.results:
        if (s & result.bitstring).bit_count() & 1:
            expectation_value -= result.count
        else:
            expectation_value += result.count
        total_count += result.count
    return expectation_value / total_count
