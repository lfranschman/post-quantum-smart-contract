import numpy as np

def mod_mult(m1, m2, q):
    return np.dot(m1, m2) % q

def mod_add(m1, m2, q):
    return (m1 + m2) % q

def mod_sub(m1, m2, q):
    return (m1 - m2) % q

def keyGen_r1cs_matrices(a, s, q):
    n = len(s)

    # Define the variable sizes
    num_s = n
    num_a = n * n
    num_e = n
    num_t = n

    # Initialize the matrices
    A = np.zeros((num_t + 2 * n, num_s + num_a + num_e))
    B = np.zeros((num_t + 2 * n, num_s + num_a + num_e))
    C = np.zeros((num_t + 2 * n, num_s + num_a + num_e))

    # Matrix `a` Constraint
    for i in range(n):
        for j in range(n):
            A[i, num_s + i * n + j] = 1
            C[i, num_t + i * n + j] = 1

    # Error Vector `e` Constraint
    for i in range(n):
        A[i + n, num_s + num_a + i] = 1
        C[i + n, num_t + i] = 1

    # Computed Value `t` Constraint
    for i in range(n):
        B[i + 2 * n, i] = 1
        B[i + 2 * n, num_s + num_a + i] = -1
        C[i + 2 * n, num_t + i] = 1

    return A, B, C

# Example usage
n = 4
q = 1666
a = np.random.randint(q, size=(n, n))
s = np.random.randint(q, size=(n, 1))
A, B, C = keyGen_r1cs_matrices(a, s, q)

print("Matrix A:")
print(A)
print("\nMatrix B:")
print(B)
print("\nMatrix C:")
print(C)
