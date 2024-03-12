import fromR1CStoQAP as qap
import TrustedSetup as ts
import numpy as np
import lwe as Kg

q = 257
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
tau, tTau, tau_a, tau_t, t_a, a, t, s, e = ts.getTau()
print("tau_a: ", tau_a) # == coefficients times tau multiple times a
print("tau_t: ", tau_t) # == coefficients times tau multiple times t

H, T = ts.hxBalancing(Ua, Va, Wa)
# [A]_1 generation
A_1 = lattice_poly_hiding(Ua, tau_a)  # sum of all tau multiples times a //times coefficient of Ua
print("A_1: ",  A_1)
# [B]_2 generation
B_2 = lattice_poly_hiding(Va, tau_t) #sum of all tau multiples times t //times coefficient of Va
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
# print("proof :", proof)
print("Ua(tau) * Va(tau) - Wa(tau) == H(T(tau)):", Ua(tau) * Va(tau) - Wa(tau) == H(tau) * T(tau))
print("Ua(tau) * Va(tau) - Wa(tau) == H(T(tau)):", Ua(tau) * Va(tau) - Wa(tau) == H(tau) * T(tau))


print("Ua(tau) * Va(tau) - Wa(tau) == H(T(tau)):", Kg.mod_mult(A_1, B_2) == np.remainder(C_1, q))
print("mult A1 B2: ", Kg.mod_mult(A_1, B_2))
print("Mult C1 HT: ", np.remainder(C_1, q))
print("a: ", a)
print("s: ", s)
print("e: ", a)
print("t: ", a)

print("check1: ", Kg.mod_mult(a, s))
print("check2: ", Kg.mod_add(Kg.mod_mult(a, s), e))
print("check3: ", np.remainder(t, q))
print(proof)
print(A_1 * B_2)
print(C_1*tTau)

# def pairing_verification(A, B, C, tTau):
#     pairing_lhs = np.dot(A, B)
#     pairing_rhs = np.dot(C, tTau)
#
#     # Check if the pairing equations hold
#     assert np.array_equal(pairing_lhs, pairing_rhs), "Pairing equations do not hold"
#
#     print("Proof verification successful")





# check = pairing_verification(A_1, B_2, C_1, t)
# print(check)
