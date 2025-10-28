from srg import srg_mielke, mielke_generator
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
from qiskit.quantum_info import SparsePauliOp
import time
import json
from pathlib import Path
from symmer import PauliwordOp
import numpy as np
import pickle
from scipy.integrate import RK45

# Load H2 Hamiltonian
folder = Path('hamiltonian_data')
ham_dict = {}
for json_file in folder.glob('*.json'):
    with open(json_file, 'r') as f:
        ham_dict[json_file.stem] = json.load(f)

mol = PauliwordOp.from_dictionary((ham_dict['H2_6-31G_SINGLET_BK']['hamiltonian']))
H = mol
H_mat = H.to_sparse_matrix
N = H_mat.shape[0]
N2 = N**2

# Data storage
num_paulis = []
H_pauli_list = []

print("Starting H2 SRG iterations...")
start_time = time.time()
H0 = H_mat.toarray()  # Convert sparse matrix to dense array
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

s = RK45(derivativeOfHt, 0.0, H0.reshape(N2), 100.0, rtol = 1e-10, atol = 1e-15)
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

iteration = 1
for i in res_t:
    H_final, _ = srg_mielke(H_mat, stepsize=i, number_of_steps=2)
    eigenvals, eigenvecs = eigsh(H_final, k=4, which='SA')
    h_pauli = SparsePauliOp.from_operator(H_final.todense())

    num_terms = len(h_pauli.paulis)
    num_paulis.append(num_terms)
    H_pauli_list.append(h_pauli.to_list())

    print(f"Iteration {iteration}: {num_terms} Pauli terms")
    iteration += 1
    H_mat = H_final

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Total execution time: {elapsed_time:.4f} seconds")

# Plot number of Pauli terms
plt.figure(figsize=(10, 6))
plt.plot(range(len(res_t)), num_paulis, 'b-', linewidth=1.5)
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
