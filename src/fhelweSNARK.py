import galois
import numpy as np
import fromLWEtoR1CS as r1
import fromR1CStoQAP as qap
import sys
import time
# import timeit


order = 73
GF = galois.GF(order)
q = 8929
t = 73
d = 57
delta = q // t
# Polynomial modulus
p_q = np.poly1d([1] + ([0] * (d - 1)) + [1])


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

add = addition(p_q, q)
mul = multiplication(p_q, q)


def mod_horner(poly, x, modulus):
    result = 0
    for coefficient in poly.coefficients:
        result = (result * x + coefficient) % modulus
    return result


def setup():
    alpha = np.random.randint(0, t, 1)[0]
    sk = np.poly1d(np.random.randint(0, 2, d))
    a2 = np.poly1d(np.random.randint(0, q, d) % q)
    e = np.poly1d(np.random.normal(0, 2, d).astype(int) % q)
    pk_0 = add(-mul(a2, sk), e)
    pk_1 = a2
    pk = (pk_0, pk_1)
    u = np.poly1d(np.random.randint(0, 2, d))
    e1 = np.poly1d(np.random.normal(0, 2, d).astype(int) % q)
    e2 = np.poly1d(np.random.normal(0, 2, d).astype(int) % q)
    return alpha, sk, a2, e, pk, u, e1, e2

def prover(pk, u, e1, alpha, s):
    A, B, C = r1.LWEToR1CS_transform()

    U, V, W, Ua, Va, Wa = qap.polySum(GF(s))

    #vanishing polynomial
    T = galois.Poly([1, order - 1], field=GF)
    for i in range(2, len(A) + 1):
        T = T * galois.Poly([1, order - i], field=GF)

    H = (Ua * Va - Wa) // T

    print("check if remainder is zero: ", (Ua * Va - Wa) % T)
    # print("Ua: ", Ua)
    # print("Va: ", Va)
    # print("Wa: ", Wa)
    # print("H: ", H)
    # print("T: ", T)
    #
    # print("Uaalpha: ", Ua(alpha))
    # print("Vaalpha: ", Va(alpha))
    # print("Waalpha: ", Wa(alpha))
    # print("Halpha: ", H(alpha))
    # print("Talpha: ", T(alpha))
    left = Ua(alpha) * Va(alpha) - Wa(alpha)
    right = H(alpha) * T(alpha)
    print("left: ", left)
    print("right: ", right)
    print("res1: ", left == right)

    c1 = add(add(mul(pk[0], u), e1), mul(delta, left))
    print("c1 check:", c1)
    c2 = add(add(mul(pk[0], u), e1), mul(delta, right))
    print("c2 check:", c2)
    # e_enc = add(mul(pk[1], u), e2)

    return c1, c2


def verifier(proof):
    c_0, c_1 = proof
    check = add(c_0, -c_1)
    print("just a check: ", int(check.coefficients[0]))

    return int(check.coefficients[0]) == 0

def main():
    s = np.array([1, 2, 4, 1, 2, 3, 2, 1, 0, 1, 1, 1, 0, 3, 0, 1, 3])
    alpha, sk, a2, e, pk, u, e1, e2 = setup()
    print(alpha)
    proof = prover(pk, u, e1, alpha, s)
    start = time.time()
    verification = verifier(proof)

    end = time.time()
    # execution_time = timeit.timeit('verifier()', globals=globals(), number=100)
    # print("average execution time over 100 runs: ", execution_time)
    print("verification result:", verification)
    print("verification time: ", end - start)
    print("proof size: ", sys.getsizeof(proof), " bytes")


if __name__ == "__main__":
    main()
    # execution_time = timeit.timeit('main()', globals=globals(), number=100)
    # print("average execution time over 100 runs: ", execution_time)
