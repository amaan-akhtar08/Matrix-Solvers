import numpy as np
import math

filename = 'input_mat/Kmat.txt'
with open(filename, 'r') as file:
    Kmat = np.loadtxt(file)

nn = int(math.sqrt(len(Kmat)))
    
A = np.reshape(Kmat, (nn, nn))

eigenvalues = np.linalg.eigvals(A)
pos_eigenvalues = abs(eigenvalues)

lambda_max = np.max(pos_eigenvalues)
lambda_min = np.min(pos_eigenvalues)
kappa = lambda_max / lambda_min

convergence_rate = (kappa - 1) / (kappa + 1)

# print(f"Eigenvalues of A: {eigenvalues}")
# print(f"Condition number (kappa): {kappa}")
print(f"#Convergence rate: {convergence_rate}")
