from srg import generate_random_pauli_sum
from srg import srg_mielke
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
from qiskit.quantum_info import SparsePauliOp
import time
import json
from pathlib import Path
from symmer import PauliwordOp
import numpy as np
folder = Path('hamiltonian_data')
ham_dict = {}
for json_file in folder.glob('*.json'):
    with open(json_file, 'r') as f:
        ham_dict[json_file.stem] = json.load(f)  # .stem gives filename without extension
mol = PauliwordOp.from_dictionary((ham_dict['H2_6-31G_SINGLET_BK']['hamiltonian']))
H = mol
H_mat = H.to_sparse_matrix
num_paulis = []
for i in range(2):
    start_time = time.time()
    H_final,_ = srg_mielke(H_mat,stepsize=0.02,number_of_steps=2)
    eigenvals,eigenvecs = eigsh(H_final,k=4,which='SA')
    h_pauli= SparsePauliOp.from_operator(H_final.todense())
    num_paulis.append(len(h_pauli.paulis))
    print(f"Number of Pauli terms after iteration {i}: {len(h_pauli.paulis)}")
    plt.scatter(i,np.abs(eigenvecs[:,0][-1]),color='blue')
    H_mat = H_final
    
    print(f"this is the {i}th iteration")
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Execution time: {elapsed_time:.4f} seconds")
