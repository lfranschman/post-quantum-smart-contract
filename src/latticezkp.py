import numpy as np
import hashlib

def hash_to_scalar(*args):
    data = ''.join(map(str, args)).encode('utf-8')
    return int(hashlib.sha256(data).hexdigest(), 16)

def generate_key_pair():
    n = 3
    q = 257
    a_first_column = np.random.randint(q, size=(n, 1))
    A = a_first_column
    for i in range(1, n):
        A = np.hstack((A, np.roll(a_first_column, i)))
    s = np.random.randint(3, size=(n, 1)) - 1
    e = np.random.randint(3, size=(n, 1)) - 1
    t = (A @ s + e) % q

    return A, t, (s, e)

def hash_challenge(*args):
    return hash_to_scalar(*args)

def prover(A, s, e):
    n = len(s)
    z = np.random.randint(0, 257, dtype=np.uint64)
    x = (A @ s + e) % 2
    y = (z * A + e) % 2  # Fix here: use A instead of s in y calculation
    c = hash_challenge(A, x, y)

    r = (z - c * s) % 257

    return c, r

def verifier(A, t, c, r):
    x = (t - c * A) % 2
    y = (r * A + c * t) % 2
    c_prime = hash_challenge(A, x, y)

    return c == c_prime

# Example usage
A, t, (s, e) = generate_key_pair()
c, r = prover(A, s, e)
result = verifier(A, t, c, r)

print("Verification Result:", result)
