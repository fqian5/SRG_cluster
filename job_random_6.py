from srg import generate_random_pauli_sum
from srg import srg_mielke
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
from qiskit.quantum_info import SparsePauliOp
import time
import numpy as np
from symmer import PauliwordOp
import json
import pickle

# Parameters
job_id = 6
seed = 1000 + job_id
n_qubits = 8
n_terms = 200
coeff_scale = 0.05

print(f"Job {job_id}: Generating random Hamiltonian with seed {seed}")
random_ham = generate_random_pauli_sum(n_qubits=n_qubits, n_terms=n_terms, coeff_scale=coeff_scale, seed=seed)
H = PauliwordOp.from_dictionary(random_ham)

H_mat = H.to_sparse_matrix

# Data storage
num_paulis = []
H_pauli_list = []

print("Starting SRG iterations...")
start_time = time.time()

for i in range(2000):
    H_final, _ = srg_mielke(H_mat, stepsize=0.02, number_of_steps=2)
    eigenvals, eigenvecs = eigsh(H_final, k=4, which='SA')
    h_pauli = SparsePauliOp.from_operator(H_final.todense())

    num_terms_current = len(h_pauli.paulis)
    num_paulis.append(num_terms_current)
    H_pauli_list.append(h_pauli.to_list())

    print(f"Iteration {i}: {num_terms_current} Pauli terms")

    H_mat = H_final

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Total execution time: {elapsed_time:.4f} seconds")

# Plot number of Pauli terms
plt.figure(figsize=(10, 6))
plt.plot(range(2000), num_paulis, 'b-', linewidth=1.5)
plt.xlabel('Iteration', fontsize=12)
plt.ylabel('Number of Pauli Terms', fontsize=12)
plt.title(f'Random Hamiltonian {job_id} (seed={seed}): Number of Pauli Terms vs SRG Iteration', fontsize=14)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f'plot_random_{job_id}.png', dpi=300)
plt.close()
print(f"Plot saved to plot_random_{job_id}.png")

# Save results
results = {
    'job_id': job_id,
    'seed': seed,
    'n_qubits': n_qubits,
    'n_terms_initial': n_terms,
    'coeff_scale': coeff_scale,
    'num_paulis': num_paulis,
    'H_pauli_list': H_pauli_list,
    'execution_time': elapsed_time
}

with open(f'results_random_{job_id}.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"Results saved to results_random_{job_id}.json with Pauli format")
