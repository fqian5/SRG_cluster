from srg import generate_random_pauli_sum
from srg import srg_mielke, mielke_generator
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
from qiskit.quantum_info import SparsePauliOp
import time
import numpy as np
from symmer import PauliwordOp
import json
import pickle
from scipy.integrate import RK45

# Parameters
job_id = 5
seed = 1000 + job_id
n_qubits = 8
n_terms = 200
coeff_scale = 0.05

print(f"Job {job_id}: Generating random Hamiltonian with seed {seed}")
random_ham = generate_random_pauli_sum(n_qubits=n_qubits, n_terms=n_terms, coeff_scale=coeff_scale, seed=seed)
H = PauliwordOp.from_dictionary(random_ham)

H_mat = H.to_sparse_matrix
N = H_mat.shape[0]
N2 = N**2

# Data storage
num_paulis = []
H_pauli_list = []

print("Starting SRG iterations...")
start_time = time.time()
H0 = H_mat.toarray()  # Convert sparse matrix to dense array
H0 = 0.5 * (H0 + H0.conj().T)  # Ensure Hermiticity
gs_support = []

def derivativeOfHt(t, H_coords):
    h = H_coords.reshape(N, N)

    # Convert to sparse matrix for mielke_generator
    import scipy.sparse as sp
    h_sparse = sp.csr_matrix(h)

    # g = np.diag(np.diag(h)) @ h - h @ np.diag(np.diag(h))
    g = mielke_generator(h_sparse)
    dh = g @ h_sparse - h_sparse @ g
    return dh.toarray().reshape(N2)

s = RK45(derivativeOfHt, 0.0, H0.reshape(N2), 1000.0, rtol = 1e-10, atol = 1e-15)
res_t = []
res_y = []
for i in range(5_000):
    res_t.append(s.t)
    res_y.append(s.y)
    if i % 1 == 0:
        # print('t_' + str(i), '=', res_t[-1])
        pass
    s.step()
    if s.status == 'finished':
        res_t.append(s.t)
        res_y.append(s.y)
        print('Finished')
        break
    if s.status == 'failed':
        print('Failed at step', i, 't =', s.t)
        break
    if i % 500 == 0:
        # e_gs_W, psi_gs_W = exact_gs_energy(s.y.reshape(N, N).real)
        # psi_gs_W = psi_gs_W.cleanup(zero_threshold = 1e-9)
        # # psi_gs_W.cleanup(zero_threshold=1e-9)
        # # psi_gs_W.sort(),
        # gs_support.append(psi_gs_W.n_terms)
        pass

# Calculate time differences
time_diffs = [res_t[i+1] - res_t[i] for i in range(len(res_t)-1)]

iteration = 1
for delta_t in time_diffs:
    H_final, _ = srg_mielke(H_mat, stepsize=delta_t, number_of_steps=2)
    eigenvals, eigenvecs = eigsh(H_final, k=4, which='SA')
    h_pauli = SparsePauliOp.from_operator(H_final.todense())
    h_pauli = h_pauli.chop(tol=1e-14)

    num_terms_current = len(h_pauli.paulis)
    num_paulis.append(num_terms_current)
    H_pauli_list.append(h_pauli.to_list())

    print(f"Iteration {iteration}: {num_terms_current} Pauli terms")
    iteration += 1
    H_mat = H_final

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Total execution time: {elapsed_time:.4f} seconds")

# Plot number of Pauli terms
plt.figure(figsize=(10, 6))
plt.plot(range(len(time_diffs)), num_paulis, 'b-', linewidth=1.5)
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
