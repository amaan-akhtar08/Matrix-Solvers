import numpy as np
import math

def calculate_spectral_radius_sor_from_file(filename,omega = 1.5):
    with open(filename, 'r') as file:
        Kmat = np.loadtxt(file)
    nn = int(math.sqrt(len(Kmat)))
    
    A = np.reshape(Kmat, (nn, nn))
    D = np.diag(np.diag(A))
    R = A-D            
    L = np.tril(A, -1)                 
    U = np.triu(A, 1)                  

    T_SOR = np.linalg.inv(D + omega * L).dot((1 - omega) * D - omega * U)
    T_J = -np.dot(np.linalg.inv(D), R)

    eigenvalues = np.linalg.eigvals(T_J)
    spectral_radius = max(abs(eigenvalues))

    with open("spectral_radius_TJ.txt", "w") as file:
        file.write(str(spectral_radius) + "\n")

    return spectral_radius

print("Python Script for Spectral Radius :")

spectral_radius = calculate_spectral_radius_sor_from_file('input_mat/Kmat.txt')

print("Spectral Radius of T_J :", spectral_radius)

conv = -(np.log(spectral_radius))
print("Convergance rate of T_J :", conv)
