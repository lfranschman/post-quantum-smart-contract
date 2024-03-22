import numpy as np
from numpy.linalg import norm

# Function to simulate discrete Gaussian sampling
def discrete_gaussian_sample(mean, sigma, dimension):
    return np.random.normal(mean, sigma, dimension).round().astype(int)

# Function to generate a random challenge from a challenge space C
def sample_challenge(C):
    return np.random.choice(C)

# Public matrices and parameters
A1 = np.array([[2, 3], [5, 7]])
A2 = np.array([[1, 1], [1, 1]])  # Example for A2, adjust as needed
Bext = np.array([[1, 0], [0, 1]])  # Example for Bext, adjust as needed

# Secret vectors and message
s1 = np.array([1, 2])
s2 = np.array([3, 4])
m = 3 # The message m

# Prover's side: Sampling y1, y2
y1 = discrete_gaussian_sample(0, 1, 2)
y2 = discrete_gaussian_sample(0, 1, 2)

m_vector = np.array([3, 0])
# Compute w
w = np.dot(A1, y1) + np.dot(A2, y2)

# Verifier's challenge
C = np.array([1, -1])  # Challenge space
c = sample_challenge(C)

# Prover's response to challenge
z1 = y1 + c * s1
z2 = y2 + c * s2

# Prover constructs t_A and t_B
t_A = A1.dot(s1) + A2.dot(s2)
t_B = Bext.dot(s2) + m

# Verifier's verification step
if norm(z1) <= 10 and norm(z2) <= 10:
    # Compute w_prime
    w_prime = A1.dot(z1) + A2.dot(z2) - c * t_A  # Adjusted t_B for challenge

    print(w)
    print(w_prime)
    # Verification check
    verification_passed = np.array_equal(w, w_prime)
else:
    verification_passed = False

print("Verification passed:", verification_passed)
