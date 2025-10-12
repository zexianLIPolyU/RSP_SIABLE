# siable
The Matlab (QCLAB) implementation of recursive state preparation (RSP) and single ancilla block encoding protocol (SIABLE) using as less CNOT gates as possible. 


## Comparison of the number of C-NOT gates between the single ancilla block encoding protocol (SIABLE) for general $2^{n-1}\times 2^{n-1}$ full-rank matrix and other unitary synthesis protocol and bounds in an $n$-qubit system.

| Number of qubits | Script | 3 | 4 | 5 | 6 | 7 | n |
| --- | --- | --- | --- | --- | --- | --- | --- |
| QSD | [qiskit_unitary_synthesis](https://quantum.cloud.ibm.com/docs/en/api/qiskit/qiskit.transpiler.passes.UnitarySynthesis) | 20 | 100 | 444 | 1868 | 7660 | $(23/48)\times4^n - (3/2)\times 2^n + (4/3)$ |
| Block-ZXZ | [test_siable_CNOT](https://github.com/zexianLIPolyU/siable/blob/main/test_siable_CNOT.m) | 19 | 95 | 423 | 1783 | 7319 | $(22/48)\times4^n - (3/2)\times 2^n + (5/3)$ |
| Shende's lower bound | - | 14 | 61 | 252 | 1020 | 4091 | $\lceil (1/4)\times(4^n - 3n - 1) \rceil$ |
| **SIABLE for full-rank matrix** <br> **(Proposed method)** | [test_siable_CNOT](https://github.com/zexianLIPolyU/siable/blob/main/test_siable_CNOT.m) | **9** | **45** | **205** | **877** | **3629** | $(11/48)\times 4^n - 2^n + (7/3)$ |
| **Lower bounds** **(Proposition 1)** | - | 6 | 29 | 125 | 508 | 2043 | $\lceil (1/8)\times4^n - (3/4)\times n \rceil$
