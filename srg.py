import scipy.sparse as sp
import numpy as np
from scipy.sparse.linalg import eigsh,expm
def mielke_generator(A):
    """
    More efficient version using sparse matrix operations
    """
    # Convert to CSR for efficient operations
    A_csr = A.tocsr()
    
    # Extract triangular parts
    lower = sp.tril(A_csr, k=-1)  # Lower triangle (excluding diagonal)
    upper = sp.triu(A_csr, k=1)   # Upper triangle (excluding diagonal)
    
    # Combine: lower unchanged + upper with flipped sign + zero diagonal
    A_modified = lower - upper
    
    return A_modified
def srg_mielke(H,stepsize=0.01,number_of_steps=100):
    srg = np.eye(H.shape[0])
    srg_list = []
    for i in range(number_of_steps):
        g = mielke_generator(H)
        u = expm(-i*stepsize * g)
        H = u @ H @ u.conj().T
        srg = u.conj().T@srg
        srg_list.append(srg)
    return H,srg_list
def generate_array(initial, final, length):
    """
    Generate an array of evenly spaced integers from initial to final.
    
    Args:
        initial: Starting value
        final: Ending value
        length: Number of elements in the array
    
    Returns:
        List of evenly spaced integers from initial to final
    """
    if length <= 0:
        raise ValueError("Length must be positive")
    
    if length == 1:
        return [initial]
    
    # Calculate the step size
    step = (final - initial) / (length - 1)
    
    # Generate the array
    result = [initial + i * step for i in range(length)]
    
    # Round to integers
    result = [int(round(x)) for x in result]
    
    return result


import numpy as np

def generate_random_pauli_sum(n_qubits, n_terms, coeff_scale=1.0, seed=None):
    """
    Generate a sum of random Pauli matrices with IID coefficients.
    
    Parameters:
    -----------
    n_qubits : int
        Number of qubits
    n_terms : int
        Number of non-zero Pauli terms to include
    coeff_scale : float
        Scale of the coefficient distribution
    seed : int, optional
        Random seed for reproducibility
        
    Returns:
    --------
    dict : Dictionary mapping Pauli strings to [real_coeff, imag_coeff]
    """
    if seed is not None:
        np.random.seed(seed)
    
    pauli_chars = ['I', 'X', 'Y', 'Z']
    pauli_sum = {}
    
    # Generate n_terms unique random Pauli strings
    while len(pauli_sum) < n_terms:
        # Randomly select Pauli operators for each qubit
        pauli_string = ''.join(np.random.choice(pauli_chars, size=n_qubits))
        
        # Only add if not already present (ensures uniqueness)
        if pauli_string not in pauli_sum:
            # Generate IID coefficient
            real_coeff = np.random.normal(0, coeff_scale)
            imag_coeff = 0.0  # Set to np.random.normal(0, coeff_scale) for complex
            
            pauli_sum[pauli_string] = [real_coeff, imag_coeff]
    
    return pauli_sum

