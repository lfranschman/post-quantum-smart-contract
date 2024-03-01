import fromR1CStoQAP as qap
import TrustedSetup as ts
import numpy as np
def lattice_poly_hiding(poly, powers_of_tau):
    coefs = poly.coefficients()
    print("coefs_check: ", coefs)
    print("len check 1: ", len(coefs))
    print("len check 2: ", len(powers_of_tau))

    assert len(coefs) == len(powers_of_tau), "Length mismatch"

    result = 0

    for i in range(len(coefs)):
        result += powers_of_tau[i] * int(coefs[-(i + 1)])

    return result


U, V, W, Ua, Va, Wa = qap.polySum()
print("Ua check: ", Ua)
tau, tTau, tau_a, tau_t, t_a, a, t = ts.getTau()

H = ts.hxBalancing(Ua, Va, Wa)
# [A]_1 generation
A_1 = lattice_poly_hiding(Ua, tau_a)
print("A_1: ",  A_1)
# [B]_2 generation
B_2 = lattice_poly_hiding(Va, tau_t)
print("B_1: ",  B_2)

# [C]_1 generation
## h(tau)t(tau)
HT_tau = lattice_poly_hiding(H, t_a)
print("HT_tau: ", HT_tau)
tmp_C = lattice_poly_hiding(Wa, tau_a)
print("tmp_C: ", tmp_C)
C_1 = tmp_C + HT_tau
print("C_1: ",  C_1)

proof = [A_1, B_2, C_1]
print("proof :", proof)
print(np.dot(A_1, B_2) == np.dot(C_1, t))
