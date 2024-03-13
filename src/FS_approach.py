import numpy as np


# Define the R1CS matrices A, B, and C
A = np.array([[2, 3], [4, 5]])
B = np.array([[1, 2], [3, 4]])
C = np.array([[8], [13]])

# Parameters for LWE
n = 2  # Dimension of LWE vectors (assuming 2 variables)
q = 257  # Modulus for LWE
alpha = 8  # Noise parameter

# Encode R1CS equations into LWE instances
# For simplicity, we'll encode each equation separately into LWE samples

# Encode the equation Ax · Bx
AB_equation = np.dot(A, B)  # Compute the product of matrices A and B
AB_sample = np.random.randint(0, q, size=n)  # Generate a random LWE sample
# Compute the right-hand side of the equation
C_expected = np.dot(AB_equation, AB_sample) % q

# Generate LWE samples
# Create an LWE instance using the parameters and the encoded equation
lwe_instance = LWE(n, q, alpha)
# Generate an LWE sample based on the right-hand side of the equation
lwe_sample = lwe_instance.sample(C_expected)

# Print the generated LWE sample
print("Generated LWE sample:", lwe_sample)