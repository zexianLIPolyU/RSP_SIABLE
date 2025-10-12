# siable
The Matlab (QCLAB) implementation of recursive state preparation (RSP) and single ancilla block encoding protocol (SIABLE) using as less CNOT gates as possible. 


## Comparison of the number of C-NOT gates between proposed recursive state preparation method (RSP) and other state preparation algorithms. 

| Methods | Script | 2 | 3 | 4 | 5 | 10 | 15 | Leading constant |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PB | [qclib](https://github.com/qclib/qclib) | 1 | 4 | 9 | 26 | 919 | 38813 | $23/24$ |
| Isometry | [qiskit](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.StatePreparation) | 1 | 4 | 11 | 26 | 1013 | 32752 | $23/24$ |
| LRSP | [qclib](https://github.com/qclib/qclib) | 1 | 4 | 9 | 21 | 913 | 30999 | $23/24$ |
| **RSP** **(Proposed method)** | [test_state_preparation](https://github.com/zexianLIPolyU/siable/blob/main/test_state_preparation.mlx) | 1 | **3** | **7** | **18** | **867** | **29627** | $11/12$ |
| **Lower bounds**| - | 1 | 2 | 5 | 12 | 505 | 16373 | $1/2$



## Comparison of the number of C-NOT gates between the single ancilla block encoding protocol (SIABLE) for general $2^{n-1}\times 2^{n-1}$ full-rank matrix and other unitary synthesis protocol and bounds in an $n$-qubit system.

| Number of qubits | Script | 3 | 4 | 5 | 6 | 7 | n |
| --- | --- | --- | --- | --- | --- | --- | --- |
| QSD | [qiskit](https://quantum.cloud.ibm.com/docs/en/api/qiskit/qiskit.transpiler.passes.UnitarySynthesis) | 20 | 100 | 444 | 1868 | 7660 | $(23/48)\times4^n - (3/2)\times 2^n + (4/3)$ |
| Block-ZXZ | [test_siable_CNOT](https://github.com/zexianLIPolyU/siable/blob/main/test_siable_CNOT.m) | 19 | 95 | 423 | 1783 | 7319 | $(22/48)\times4^n - (3/2)\times 2^n + (5/3)$ |
| Shende's lower bound | - | 14 | 61 | 252 | 1020 | 4091 | $\lceil (1/4)\times(4^n - 3n - 1) \rceil$ |
| **SIABLE for full-rank matrix** <br> **(Proposed method)** | [test_siable_CNOT](https://github.com/zexianLIPolyU/siable/blob/main/test_siable_CNOT.m) | **9** | **45** | **205** | **877** | **3629** | $(11/48)\times 4^n - 2^n + (7/3)$ |
| **Lower bounds**| - | 6 | 29 | 125 | 508 | 2043 | $\lceil (1/8)\times4^n - (3/4)\times n \rceil$
