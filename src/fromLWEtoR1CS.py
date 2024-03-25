import numpy as np
import sys
import lwe
import fromR1CStoQAP as qap
import galois
from numpy.testing import assert_array_equal

order = 257
GF = galois.GF(order)

np.set_printoptions(threshold=sys.maxsize)

n = 2
t = 17
# Highest coefficient power used
d = 8
# Coefficient modulus
q = 877
delta = q // t
# Polynomial modulus
p_q = np.poly1d([1] + ([0] * (d - 1)) + [1])
# a * s + e = t


def LWEToR1CS_transform():

    num_gates = n * n * 2
    num_variables = n * n + 2 * n + num_gates + 1

    # Initialize matrices with zeros
    A_matrix = np.zeros((num_gates, num_variables), dtype=int)
    B_matrix = np.zeros((num_gates, num_variables), dtype=int)
    C_matrix = np.zeros((num_gates, num_variables), dtype=int)

    #format will be: {1, t0-tn (n), A00 - Ann(n^2), s0 - sn(n), e0 - en(n), gates a0 - an^2 (n^2), gates c0-c n^2 - n (n^2 - n)}
    #processing gates a0 - an^2
    for i in range(n*n):
        A_matrix[i][i + n + 1] = 1
        B_matrix[i][1 + n + n*n + (i % n)] = 1
        C_matrix[i][1 + 3*n + n*n + i] = 1

    # processing gates c0-cn^2 and outputs/t
    c = 0
    d = 0
    z = 0
    p = 0
    for i in range(n*n):
        B_matrix[n*n + i][0] = 1
        if (i % n) == 0:
            A_matrix[n*n + i][1 + 3*n + i + n*n] = 1
            A_matrix[n*n + i][1 + 3*n + n*n + i + 1] = 1
            C_matrix[n*n + i][1 + 3 * n + z + 2 * n * n] = 1
            z = z + 1

        elif ((i + 1) % n) != 0:
            A_matrix[n * n + i][1 + 3 * n + i + n * n + 1] = 1
            A_matrix[n*n + i][1 + 3 * n + c + 2 * n * n] = 1
            C_matrix[n*n + i][1 + 3 * n + z + 2 * n * n] = 1
            c = c + 1
            z = z + 1

        else:
            if 1 + 3 * n + c + 2 * n * n < num_variables:
                A_matrix[n*n + i][1 + 3 * n + c + 2 * n * n] = 1
                A_matrix[n*n + i][1 + 2 * n + n*n + d] = 1
                C_matrix[n*n + i][1 + p] = 1
                c = c + 1
                d = d + 1
                p = p + 1

    return A_matrix, B_matrix, C_matrix


def polymod(poly, poly_mod, coeff_mod):
    """
    Computes the remainder after a polynomial division
    Args:
        poly: Polynomial
        poly_mod: Polynomial modulus
        coeff_mod: Coefficient modulus
    Returns:
        The coefficients of the remainder when `poly` is divided by `poly_mod`
    """
    return np.poly1d(np.floor(np.polydiv(poly, poly_mod)[1]) % coeff_mod)


def addition(poly_mod, coeff_mod):
    """
    Creates a function which performs polynomial addition and auto-applys polynomial- and coefficient modulus
    Args:
        poly_mod: Polynomial modulus
        coeff_mod: Coefficient modulus
    Returns:
        A function which takes polynomials `a` and `b` and adds them together
    """
    return lambda a, b: np.poly1d(polymod(np.polyadd(a, b), poly_mod, coeff_mod))


def multiplication(poly_mod, coeff_mod):
    """
    Creates a function which performs polynomial multiplication and auto-applys polynomial- and coefficient modulus
    Args:
        poly_mod: Polynomial modulus
        coeff_mod: Coefficient modulus
    Returns:
        A function which takes polynomials `a` and `b` and multiplies them
    """
    return lambda a, b: np.poly1d(polymod(np.polymul(a, b), poly_mod, coeff_mod))


def setup():
    e = np.poly1d(np.random.normal(0, 2, d).astype(int) % q)
    a = np.poly1d(np.random.randint(0, q, d) % q)
    vk = np.poly1d(np.random.randint(0, 2, d))
    alpha = np.random.randint(t)
    return e, a, vk, alpha

def mod_horner(poly, x, modulus):
    result = 0
    for coefficient in poly.coefficients:
        result = (result * x + coefficient) % modulus
    return result

def construct_proof(alpha):
    # Evaluate A, B, C polynomials at alpha
    print("alpha: ", alpha)
    A, B, C, As, Bs, Cs = qap.polySum()
    # print("U: ", A)
    # print("V: ", B)
    # print("W: ", C)
    #
    # print("Ua: ", As)
    # print("Va: ", Bs)
    # print("Wa: ", Cs)
    T_coefficients = [1]
    T = galois.Poly(T_coefficients, field=galois.GF(order))
    for i in range(2, len(A) + 1):
        T *= galois.Poly([1, -i], field=galois.GF(order))
    # Construct the solution polynomial H(x) from the witness and evaluate it at alpha
    H = (As * Bs - Cs) // T
    # H_alpha = H(alpha)
    left_side = As(alpha) * Bs(alpha) - Cs(alpha)
    right_side = H(alpha) * T(alpha)
    print(As)
    print(Bs)
    print(Cs)
    print(H)
    print(T)
    print(left_side)
    print(right_side)
    return left_side == right_side

    # add = addition(p_q, q)
    # mul = multiplication(p_q, q)
    # #
    # u = np.poly1d(np.random.randint(0, 2, d))
    # e_1 = np.poly1d(np.random.normal(0, 1, d).astype(int) % q)
    # e_2 = np.poly1d(np.random.normal(0, 1, d).astype(int) % q)

    # print("pk: ", pk)
    # print("u: ", u)
    # print("e_1: ", e_1)
    # print("delta: ", delta)

    # Ascoefs = As.coefficients()
    # As1d = np.poly1d(np.array(Ascoefs))
    # print("As: ", As)
    # print("As1d: ", As1d)
    # print("As_alpha: ", As(alpha))
    # print("As1d_alpha: ", mod_horner(As1d, alpha, q))
    #
    # Bscoefs = Bs.coefficients()
    # Bs1d = np.poly1d(np.array(Bscoefs))
    # print("Bs: ", Bs)
    # print("Bs1d: ", Bs1d)
    # print("Bs_alpha: ", Bs(alpha))
    # print("Bs1d_alpha: ", mod_horner(Bs1d, alpha, q))
    #
    # Cscoefs = Cs.coefficients()
    # Cs1d = np.poly1d(np.array(Cscoefs))
    # print("Cs: ", Cs)
    # print("Cs1d: ", Cs1d)
    # print("Cs_alpha: ", Cs(alpha))
    # print("Cs1d_alpha: ", mod_horner(Cs1d, alpha, q))
    #
    # Tcoefs = T.coefficients()
    # T1d = np.poly1d(np.array(Tcoefs))
    # print("T: ", T)
    # print("T1d: ", T1d)
    # print("T_alpha: ", T(alpha))
    # print("T1d_alpha: ", mod_horner(T1d, alpha, q))
    #
    # Hcoefs = H.coefficients()
    # H1d = np.poly1d(np.array(Hcoefs))
    # print("H: ", H)
    # print("H1d: ", H1d)
    # print("H_alpha: ", H(alpha))
    # print("H1d_alpha: ", mod_horner(H1d, alpha, q))
    #
    # m = np.poly1d((np.array([4, 3, 6, 7, 7, 8, 3, 2])))
    # print("m: ", m)
    # # # comb_AB = mul(As1d, Bs1d)
    # e_enc = add(mul(pk[1], u), e_2)
    # A_enc = add(add(mul(pk[0], u), e_1), mul(delta, m))
    # B_enc = add(add(mul(pk[0], u), e_1), mul(delta, Bs1d))
    # C_enc = add(add(mul(pk[0], u), e_1), mul(delta, Cs1d))
    # T_enc = add(add(mul(pk[0], u), e_1), mul(delta, T1d))
    # H_enc = add(add(mul(pk[0], u), e_1), mul(delta, H1d))

    # A_prime = np.poly1d(np.round(add(mul(e_enc, vk), A_enc) * t / q))
    # print("primeCheck: ", A_prime)
    # print("A_primeeeeee: ", A_prime)
    #
    # As1d_alpha = mod_horner(As1d, alpha, q)
    # Bs1d_alpha = mod_horner(Bs1d, alpha, q)
    # Cs1d_alpha = mod_horner(Cs1d, alpha, q)
    # T1d_alpha = mod_horner(T1d, alpha, q)
    # H1d_alpha =  mod_horner(H1d, alpha, q)



    # enc_As1d_alpha = mod_horner(A_enc, alpha, q)
    # enc_Bs1d_alpha = mod_horner(B_enc, alpha, q)
    # enc_Cs1d_alpha = mod_horner(C_enc, alpha, q)
    # enc_T1d_alpha = mod_horner(T_enc, alpha, q)
    # enc_H1d_alpha =  mod_horner(H_enc, alpha, q)


    # left_side = As_alpha * Bs_alpha - Cs_alpha
    # right_side = H_alpha * T(alpha)
    # assert left_side == right_side, f"QAP relation does not hold: {left_side} != {right_side}"
    # print("it passed!!")
    #
    # left_side = (As1d_alpha * Bs1d_alpha - Cs1d_alpha) % q
    # right_side = (T1d_alpha * H1d_alpha) % q
    # assert left_side == right_side, f"QAP relation does not hold: {left_side} != {right_side}"
    # print("it passed!!")
    #
    # # left_side_encrypted = enc_As1d_alpha * enc_Bs1d_alpha - enc_Cs1d_alpha
    # # right_side_encrypted = enc_T1d_alpha * enc_H1d_alpha
    # # assert abs(left_side_encrypted - right_side_encrypted) < .1 , f"QAP relation does not hold: {left_side_encrypted} != {right_side_encrypted}"
    # # print("it passed!!")
    # #
    # # ab_comb = As(alpha) * Bs(alpha)
    # # print(ab_comb - Cs(alpha))
    # # print(T(alpha))
    # # print(H(alpha))
    # #
    # # print(public_key)
    # # e_2 = np.poly1d(np.random.normal(0, 2, d).astype(int) % q)
    #
    # return {'enc_e': e_enc, 'enc_H': H_enc, 'enc_T': T_enc,  'enc_A': A_enc, 'enc_B': B_enc, 'enc_C': C_enc}





def verify_proof(proof, gate_checks, vk, s, a, alpha, public_key):
    # if np.all(gate_checks) == False:
    #     return False
    # Extract the encrypted evaluations from the proof
    enc_e = proof['enc_e']
    enc_A = proof['enc_A']
    enc_B = proof['enc_B']
    enc_C = proof['enc_C']
    enc_T = proof['enc_T']
    enc_H = proof['enc_H']

    add = addition(p_q, q)
    mul = multiplication(p_q, q)

    # u1, v1 = enc_left
    # u2, v2 = enc_T
    # u3, v3 = enc_H
    # A_prime = np.poly1d(np.round(add(mul(enc_e, vk), enc_A) * t / q) % t)
    # print("A_prime: ", A_prime)
    A = np.poly1d(np.round(add(mul(enc_e, vk), enc_A) * t / q) % t)
    #lwe.decrypt(s, u1,  v1)
    B =  np.poly1d(np.round(add(mul(enc_e, vk), enc_B) * t / q) % t)
    #lwe.decrypt(s, u2, v2)
    C = np.poly1d(np.round(add(mul(enc_e, vk), enc_C) * t / q) % t)
    T = np.poly1d(np.round(add(mul(enc_e, vk), enc_T) * t / q) % t)
    H = np.poly1d(np.round(add(mul(enc_e, vk), enc_H) * t / q) % t)
    #lwe.decrypt(s, u3, v2)
    print("A: ", A)
    print("B: ", B)
    print("C: ", C)
    print("T: ", T)
    print("H: ", H)

    # left_side = p1
    # right_side = p2 * p3
    # assert left_side == right_side, f"QAP relation does not hold: {left_side} != {right_side}"
    # print("it passed!!")
    # print(enc_A_alpha)
    # print(enc_B_alpha)
    # print(enc_C_alpha)
    # left_side_encrypted = enc_AB_alpha - enc_C_alpha
    # right_side_encrypted = enc_T_alpha * enc_H_alpha
    # # print(enc_left_side)
    # # print(enc_right_side)
    # qap_check = left_side_encrypted - right_side_encrypted < .1

    add = addition(p_q, q)
    mul = multiplication(p_q, q)
    # The proof is valid if the encrypted left side is close enough to the encrypted right side, considering the encryption noise
    return A(alpha) * B(alpha) - C(alpha) == T(alpha) * H(alpha)


def main():

    #get the r1cs
    A, B, C = LWEToR1CS_transform()
    assert q == delta * t + (q % t), "does not hold"
    assert p_q.order == d, "does not hold"
    witness = np.array([1, 22, 38, 21, 14, 37, 19, 1, 0, 1, 1, 21, 0, 37, 0, 21, 37])
    # print(witness)
    print(np.matmul(C, witness) == (np.matmul(A, witness) * np.matmul(B, witness)))
    gate_checks = np.matmul(C, witness) == (np.matmul(A, witness) * np.matmul(B, witness))
    #generate lwe public key
    # a_first_column = np.random.randint(q, size=(n, 1))
    # m = 2
    # a = a_first_column
    # for i in range(1, m):
    #     a = np.hstack((a, np.roll(a_first_column, i)))
    # s = np.random.randint(3, size=(m, 1)) - 1
    # e = np.random.randint(3, size=(m, 1)) - 1
    # public_key = lwe.keyGen(a, s, e)
    # print(public_key)
    #get the QAP from the R1CS
    A_polys = qap.get_polys_of_matrix(A)
    B_polys = qap.get_polys_of_matrix(B)
    C_polys = qap.get_polys_of_matrix(C)

    # #generate the target polynomial
    # T_coefficients = [1]
    # T = galois.Poly(T_coefficients, field=galois.GF(order))
    # for i in range(2, len(A) + 1):
    #     T *= galois.Poly([1, -i], field=galois.GF(order))
    # print("witness: ", witness)




    e, a, vk, alpha = setup()
    add = addition(p_q, q)
    mul = multiplication(p_q, q)
    pk0 = add(-mul(a, vk), e)
    pk1 = a
    pk = (pk0, pk1)
    extr_e = add(mul(pk[1], vk), pk[0])
    assert_array_equal(extr_e, e)
    proof = construct_proof(alpha)
    print(proof)
    # new_gate_checks = np.append(gate_checks, np.append(gate_checks, gate_checks))
    # verifiction = verify_proof(proof, new_gate_checks, vk, a, alpha)
    # print(verifiction)



if __name__ == "__main__":
    main()
