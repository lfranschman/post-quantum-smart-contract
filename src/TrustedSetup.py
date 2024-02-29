import galois
import numpy as np
import fromLWEtoR1CS as r1
import KeyGen as Kg
import fromR1CStoQAP as qap

order = 2**8 + 1
GF = galois.GF(order)
n = 2
q = 257

a_first_column = np.random.randint(q, size=(n, 1))
a = a_first_column
for i in range(1, n):
    a = np.hstack((a, np.roll(a_first_column, i)))
s = np.random.randint(3, size=(n, 1)) - 1
e = np.random.randint(3, size=(n, 1)) - 1

t = Kg.keyGen(a, s, e)


# tau value is random
tau = np.random.randint(q)
# alpha value is random
alpha = np.random.randint(q)
# beta value is random
beta = np.random.randint(q)


# n degree of polynomials
A, B, C = r1.LWEToR1CS_transform()
n = len(A)-1
# public inputs (matrix a and vector t)
l = 2
# private inputs = total inputs - public inputs
m = len(A[0])

# Compute T(x)
T_coefficients = [-1] + [0] * (n - 1)  # T = (x - 1)
T = galois.Poly(T_coefficients, field=galois.GF(q))
for i in range(2, n + 1):
    T *= galois.Poly([-i, 1], field=galois.GF(q))

# t(tau)
tTau = T(tau)

# Powers of tau for [A] and [C] (on lattice)
tau_A = A * tau
tau_t = t * tau

# Powers of tau for h(tau)t(tau)
t_A = A * tTau
t_t = t * tTau

# Print the results
print("tau: ", tau)
print("Ttau: ", tTau)
print("Powers of tau for [A]:\n", tau_A)
print("Powers of tau for [t]:\n", tau_t)
print("Powers of tau for h(tau)t(tau) [A]:\n", t_A)
print("Powers of tau for h(tau)t(tau) [t]:\n", t_t)


def hxBalancing(Ua, Va, Wa):
    H = (Ua * Va - Wa) // T
    print("t(x):", T)
    print("h(x):", H)
    # remainder should be null
    assert (Ua * Va - Wa) % T == 0, "Remainder not zero!"
    print("Ua(tau) * Va(tau) - Wa(tau) == H(T(tau)):", Ua(tau) * Va(tau) - Wa(tau) == H(tau) * T(tau))
    print("deg(T) =", T.degree)
    print("deg(H) =", H.degree)


U, V, W, Ua, Va, Wa = qap.polySum()
hxBalancing(Ua, Va, Wa)


