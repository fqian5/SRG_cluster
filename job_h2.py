from srg import srg_mielke
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
from qiskit.quantum_info import SparsePauliOp
import time
import json
from pathlib import Path
from symmer import PauliwordOp
import numpy as np
import pickle

# Load H2 Hamiltonian
folder = Path('hamiltonian_data')
ham_dict = {}
for json_file in folder.glob('*.json'):
    with open(json_file, 'r') as f:
        ham_dict[json_file.stem] = json.load(f)

mol = PauliwordOp.from_dictionary((ham_dict['H2_6-31G_SINGLET_BK']['hamiltonian']))
H = mol
H_mat = H.to_sparse_matrix

# Data storage
num_paulis = []
H_pauli_list = []

print("Starting H2 SRG iterations...")
start_time = time.time()

for i in range(2000):
    H_final, _ = srg_mielke(H_mat, stepsize=0.02, number_of_steps=2)
    eigenvals, eigenvecs = eigsh(H_final, k=4, which='SA')
    h_pauli = SparsePauliOp.from_operator(H_final.todense())

    num_terms = len(h_pauli.paulis)
    num_paulis.append(num_terms)
    H_pauli_list.append(h_pauli.to_list())

    print(f"Iteration {i}: {num_terms} Pauli terms")

    H_mat = H_final

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Total execution time: {elapsed_time:.4f} seconds")

# Plot number of Pauli terms
plt.figure(figsize=(10, 6))
plt.plot(range(2000), num_paulis, 'b-', linewidth=1.5)
plt.xlabel('Iteration', fontsize=12)
plt.ylabel('Number of Pauli Terms', fontsize=12)
plt.title('H2: Number of Pauli Terms vs SRG Iteration', fontsize=14)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('plot_h2.png', dpi=300)
plt.close()
print("Plot saved to plot_h2.png")

# Save results
results = {
    'molecule': 'H2_6-31G_SINGLET_BK',
    'num_paulis': num_paulis,
    'H_pauli_list': H_pauli_list,
    'execution_time': elapsed_time
}

with open('results_h2.json', 'w') as f:
    json.dump(results, f, indent=2)

print("Results saved to results_h2.json with Pauli format")
