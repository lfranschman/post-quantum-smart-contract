import numpy as np
import sys
import lwe
import fromR1CStoQAP as qap
import galois

order = 887
GF = galois.GF(order)

np.set_printoptions(threshold=sys.maxsize)

n = 2
q = 887
d = 2 ** n
t = 7
delta = q // t
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
    sk = np.poly1d(np.random.randint(0, 2, d))
    alpha = np.random.randint(q)
    return e, a, sk, alpha


def construct_proof(alpha, pk, witness, T):
    # Evaluate A, B, C polynomials at alpha
    A, B, C, As, Bs, Cs = qap.polySum(witness)
    # print("U: ", A)
    # print("V: ", B)
    # print("W: ", C)
    #
    # print("Ua: ", As)
    # print("Va: ", Bs)
    # print("Wa: ", Cs)
    As_alpha = As(alpha)
    Bs_alpha = Bs(alpha)
    Cs_alpha = Cs(alpha)


    # Construct the solution polynomial H(x) from the witness and evaluate it at alpha
    H = (As * Bs - Cs) // T
    H_alpha = H(alpha)

    add = addition(p_q, q)
    mul = multiplication(p_q, q)

    u = np.poly1d(np.random.randint(0, 2, d))
    e_1 = np.poly1d(np.random.normal(0, 1, d).astype(int) % q)
    # print("pk: ", pk)
    # print("u: ", u)
    # print("e_1: ", e_1)
    # print("delta: ", delta)
    Ascoefs = As.coefficients()
    As1d = np.poly1d(Ascoefs)
    Bscoefs = Bs.coefficients()
    Bs1d = np.poly1d(Bscoefs)
    Cscoefs = Cs.coefficients()
    Cs1d = np.poly1d(Cscoefs)
    Tcoefs = T.coefficients()
    T1d = np.poly1d(Tcoefs)
    Hcoefs = H.coefficients()
    H1d = np.poly1d(Hcoefs)

    comb_AB = mul(As1d, Bs1d)
    # c_0 = add(add(mul(pk, u), e_1), mul(delta, As1d))
    c_1 = add(add(mul(pk, u), e_1), mul(delta, comb_AB))
    c_2 = add(add(mul(pk, u), e_1), mul(delta, Cs1d))
    c_3 = add(add(mul(pk, u), e_1), mul(delta, T1d))
    c_4 = add(add(mul(pk, u), e_1), mul(delta, H1d))

    # print(c_0)


    # # QAP relation check at alpha
    # left_side = As_alpha * Bs_alpha - Cs_alpha
    # right_side = H_alpha * T(alpha)
    # assert left_side == right_side, f"QAP relation does not hold: {left_side} != {right_side}"
    # print("it passed!!")
    #
    # left_side_encrypted = c_1(alpha) - c_2(alpha)
    # right_side_encrypted = c_3(alpha) * c_4(alpha)
    # assert (left_side_encrypted - right_side_encrypted) < 1.1 , f"QAP relation does not hold: {left_side_encrypted} != {right_side_encrypted}"
    # print("it passed!!")

    return {'enc_H_alpha': c_4(alpha), 'enc_T_alpha': c_3(alpha),  'enc_AB_alpha': c_1(alpha), 'enc_C_alpha': c_2(alpha)}





def verify_proof(proof):
    # Extract the encrypted evaluations from the proof
    enc_AB_alpha = proof['enc_AB_alpha']
    enc_C_alpha = proof['enc_C_alpha']
    enc_T_alpha = proof['enc_T_alpha']
    enc_H_alpha = proof['enc_H_alpha']

    # print(enc_A_alpha)
    # print(enc_B_alpha)
    # print(enc_C_alpha)
    left_side_encrypted = enc_AB_alpha - enc_C_alpha
    right_side_encrypted = enc_T_alpha * enc_H_alpha
    # print(enc_left_side)
    # print(enc_right_side)

    # The proof is valid if the encrypted left side is close enough to the encrypted right side, considering the encryption noise
    # Adjust the comparison as per your encryption scheme's noise characteristics
    return left_side_encrypted - right_side_encrypted < .1


def main():

    #get the r1cs
    A, B, C = LWEToR1CS_transform()
    witness = np.array([1, 22, 38, 21, 14, 37, 19, 1, 0, 1, 1, 21, 0, 37, 0, 21, 37])
    # print(witness)
    print(np.matmul(C, witness) == (np.matmul(A, witness) * np.matmul(B, witness)))

    #generate lwe public key
    a_first_column = np.random.randint(q, size=(n, 1))
    a = a_first_column
    for i in range(1, n):
        a = np.hstack((a, np.roll(a_first_column, i)))
    s = np.random.randint(3, size=(n, 1)) - 1
    e = np.random.randint(3, size=(n, 1)) - 1
    public_key = lwe.keyGen(a, s, e)

    #get the QAP from the R1CS
    A_polys = qap.get_polys_of_matrix(A)
    B_polys = qap.get_polys_of_matrix(B)
    C_polys = qap.get_polys_of_matrix(C)

    #generate the target polynomial
    T_coefficients = [1]
    T = galois.Poly(T_coefficients, field=galois.GF(q))
    for i in range(2, len(A) + 1):
        T *= galois.Poly([1, -i], field=galois.GF(q))
    # print("witness: ", witness)

    e, a, sk, alpha = setup()
    add = addition(p_q, q)
    mul = multiplication(p_q, q)
    pk = add(-mul(a, sk), e)
    proof = construct_proof(alpha, pk, GF(witness), T)
    verifiction = verify_proof(proof)
    print(verifiction)



if __name__ == "__main__":
    main()
