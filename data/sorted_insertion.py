# Everything in this file is based on: https://doi.org/10.22331/q-2021-01-20-385

from typing import Dict, List


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
        denominator += term**0.5
    return (numerator / denominator) ** 2


def are_commuting(p, q):
    """
    p, q: Pauli strings of same length
    """
    sign = 1
    for i in range(len(p)):
        if p[i] == "I" or q[i] == "I" or p[i] == q[i]:
            pass
        else:
            sign *= -1
    return sign == 1


def commutes_with_all(p, paulis):
    for q in paulis:
        if not are_commuting(p, q):
            return False
    return True


def are_qubitwise_commuting(p, q):
    """
    p, q: Pauli strings of same length
    """
    for i in range(len(p)):
        if not (p[i] == "I" or q[i] == "I" or p[i] == q[i]):
            return False
    return True


def qubitwise_commutes_with_all(p, paulis):
    for q in paulis:
        if not are_qubitwise_commuting(p, q):
            return False
    return True


def sorted_insertion_qwc(ham):
    # sort Paulis: biggest coefficient (abs value) first
    sorted_paulis = {
        k: v for k, v in sorted(ham.items(), key=lambda item: abs(item[1]))[::-1]
    }
    out = []
    out.append({"operators": []})
    i_max = 1  # number of sets
    for p in sorted_paulis:  # Loop through all Pauli operators p
        i = 0
        assigned = False
        while not assigned:  # Find a set into which p fits
            if qubitwise_commutes_with_all(p, out[i]["operators"]):
                out[i]["operators"].append(p)
                assigned = True
            else:
                i += 1
                if (
                    i == i_max
                ):  # p did not fit into any existing set --> create a new one
                    out.append({"operators": []})
                    i_max = i + 1
    return out


def sorted_insertion(ham):
    # sort Paulis: biggest coefficient (abs value) first
    sorted_paulis = {
        k: v for k, v in sorted(ham.items(), key=lambda item: abs(item[1]))[::-1]
    }
    out = []
    out.append({"operators": []})
    i_max = 1  # number of sets
    for p in sorted_paulis:  # Loop through all Pauli operators p
        i = 0
        assigned = False
        while not assigned:  # Find a set into which p fits
            if commutes_with_all(p, out[i]["operators"]):
                out[i]["operators"].append(p)
                assigned = True
            else:
                i += 1
                if (
                    i == i_max
                ):  # p did not fit into any existing set --> create a new one
                    out.append({"operators": []})
                    i_max = i + 1
    return out
