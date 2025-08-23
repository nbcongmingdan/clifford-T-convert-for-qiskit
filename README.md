# clifford-T-convert-for-qiskit
A python program to translate qiskit quantum circuit into fault tolerent style(clifford+T) circuit.

# Quantum GridSynthesis Demo

This small application depends on:

- [Qiskit](https://qiskit.org/)  
- The **PyGridsynth** algorithm by Shuntaro Yamamoto and Nobuyuki Yoshioka  
  Repository: https://github.com/quantum-programming/pygridsynth

## Installation

Install the required packages via pip:

```bash
pip install qiskit pygridsynth

```

Usage

```python

from clifford_t_converter import to_clifford_t

# qc_test is a Qiskit QuantumCircuit object
qc_conv = to_clifford_t(qc_test, precision=1e-3)

```
