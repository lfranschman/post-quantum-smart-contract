import numpy as np

# Parameters
q = 257
n = 2

def generate_keypair():
    A = np.random.randint(q, size=(n, n))
    s = np.random.randint(3, size=(n, 1)) - 1
    e = np.random.randint(3, size=(n, 1)) - 1

    t = (A.dot(s) + e) % q

    return A, t,e, s

def prove_knowledge(s):
    # Simulate a proof by simply returning the secret vector
    return s

def verify_knowledge(A, t, e, proof):
    # Verify the knowledge by checking if A * proof is equal to t
    result = (A.dot(proof) % q) + e
    return np.array_equal(result, t)

# Example usage
A, t, e, s = generate_keypair()
print("Matrix A:")
print(A)
print("Target Vector t:")
print(t)
print("Secret Vector s:")
print(s)

proof = prove_knowledge(s)
print("Proof:")
print(proof)

verification_result = verify_knowledge(A, t, e, proof)
print("Verification Result:", verification_result)
