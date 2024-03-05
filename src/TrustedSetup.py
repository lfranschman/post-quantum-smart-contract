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
print("n check: ", n)
# public inputs (matrix a and vector t)
l = 2
# private inputs = total inputs - public inputs
m = len(A[0])

T_coefficients = [1]
T = galois.Poly(T_coefficients, field=galois.GF(q))
for i in range(2, len(A) + 1):
    T *= galois.Poly([1, -i], field=galois.GF(q))

def getTau():
    # t(tau)
    tTau = T(tau)

    # powers of tau for [A] and [C] (on a)
    tau_a = []
    # powers of tau for [B] (on t)
    tau_t = []
    # Powers of tau for [A] and [C] (on lattice)
    for i in range(T.degree + 1):
        local_tau = int(tau ** i)
        tau_a.append(a * local_tau)
        tau_t.append(t * local_tau)

    # Powers of tau for h(tau)t(tau)
    t_a = []
    for i in range(T.degree + 1):
        local_tau = int(tau ** i)
        tmp = (local_tau * T(tau))
        t_a.append(a * int(tmp))
    return tau, tTau, tau_a, tau_t, t_a, a, t, s, e


tau, tTau, tau_A, tau_t, t_A, a, t, s, e = getTau()
# Print the results
# print("tau: ", tau)
# print("Ttau: ", tTau)
# print("Powers of tau for [A]:\n", tau_A)
# print("Powers of tau for [t]:\n", tau_t)
# print("Powers of tau for h(tau)t(tau) [A]:\n", t_A)
# print("Powers of tau for h(tau)t(tau) [t]:\n", t_t)


def hxBalancing(Ua, Va, Wa):
    H = (Ua * Va - Wa) // T
    print("t(x):", T)
    print("h(x):", H)
    # remainder should be null
    assert (Ua * Va - Wa) % T == 0, "Remainder not zero!"
    print("Ua(tau) * Va(tau) - Wa(tau) == H(T(tau)):", Ua(tau) * Va(tau) - Wa(tau) == H(tau) * T(tau))
    print("it is: ", Ua(tau) * Va(tau) - Wa(tau))
    print("deg(T) =", T.degree)
    print("deg(H) =", H.degree)
    return H


U, V, W, Ua, Va, Wa = qap.polySum()
hxBalancing(Ua, Va, Wa)


