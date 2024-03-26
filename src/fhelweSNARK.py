import galois
import numpy as np
import fromLWEtoR1CS as r1
import fromR1CStoQAP as qap


order = 73
GF = galois.GF(order)
q = 8929
t = 73
d = 15
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

def to_poly1d(As, Bs, Cs, H, T, alpha):
    Ascoefs = As.coefficients()
    As1d = np.poly1d(np.array(Ascoefs))
    print("As: ", As)
    print("As1d: ", As1d)
    print("As_alpha: ", As(alpha))
    print("As1d_alpha: ", mod_horner(As1d, alpha, t))

    Bscoefs = Bs.coefficients()
    Bs1d = np.poly1d(np.array(Bscoefs))
    print("Bs: ", Bs)
    print("Bs1d: ", Bs1d)
    print("Bs_alpha: ", Bs(alpha))
    print("Bs1d_alpha: ", mod_horner(Bs1d, alpha, t))

    Cscoefs = Cs.coefficients()
    Cs1d = np.poly1d(np.array(Cscoefs))
    print("Cs: ", Cs)
    print("Cs1d: ", Cs1d)
    print("Cs_alpha: ", Cs(alpha))
    print("Cs1d_alpha: ", mod_horner(Cs1d, alpha, t))

    Tcoefs = T.coefficients()
    T1d = np.poly1d(np.array(Tcoefs))
    print("T: ", T)
    print("T1d: ", T1d)
    print("T_alpha: ", T(alpha))
    print("T1d_alpha: ", mod_horner(T1d, alpha, t))

    Hcoefs = H.coefficients()
    H1d = np.poly1d(np.array(Hcoefs))
    print("H: ", H)
    print("H1d: ", H1d)
    print("H_alpha: ", H(alpha))
    print("H1d_alpha: ", mod_horner(H1d, alpha, t))
    return As1d, Bs1d, Cs1d, T1d, H1d

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

def prover(pk, u, e1, e2, alpha):
    A, B, C = r1.LWEToR1CS_transform()

    U, V, W, Ua, Va, Wa = qap.polySum()
    T = galois.Poly([1, order - 1], field=GF)
    ## Then multiply by the other values: (x-2)(x-3)...(x-7)
    for i in range(2, len(A) + 1):
        T = T * galois.Poly([1, order - i], field=GF)

    H = (Ua * Va - Wa) // T

    print("ok lets see: ", (Ua * Va - Wa) % T)
    print("Ua: ", Ua)
    print("Va: ", Va)
    print("Wa: ", Wa)
    print("H: ", H)
    print("T: ", T)

    print("Uaalpha: ", Ua(alpha))
    print("Vaalpha: ", Va(alpha))
    print("Waalpha: ", Wa(alpha))
    print("Halpha: ", H(alpha))
    print("Talpha: ", T(alpha))
    left = (Ua(alpha) * Va(alpha)) - Wa(alpha)
    right = H(alpha) * T(alpha)
    print("left: ", left)
    print("right: ", right)
    print("res1: ", left == right)

    a, b, c, t2, h = to_poly1d(Ua, Va, Wa, T, H, alpha)
    print("a: ", a)
    print("b: ", b)
    print("c: ", c)
    print("t2: ", t2)
    print("h: ", h)
    a_alpha = mod_horner(a, alpha, t)
    b_alpha = mod_horner(b, alpha, t)
    c_alpha = mod_horner(c, alpha, t)
    t_alpha = mod_horner(t2, alpha, t)
    h_alpha = mod_horner(h, alpha, t)
    print("a_alpha: ", a_alpha)
    print("b_alpha: ", b_alpha)
    print("c_alpha: ", c_alpha)
    print("t_alpha: ", t_alpha)
    print("h_alpha: ", h_alpha)
    left2 = (a_alpha * b_alpha - c_alpha) % t
    right2 = (t_alpha * h_alpha) % t
    print("left2: ", left2)
    print("right2: ", right2)
    print("res2: ", left2 == right2)

    e_enc = add(mul(pk[1], u), e2)
    c0 = add(add(mul(pk[0], u), e1), mul(delta, a))
    c1 = add(add(mul(pk[0], u), e1), mul(delta, b))
    c2 = add(add(mul(pk[0], u), e1), mul(delta, c))
    c3 = add(add(mul(pk[0], u), e1), mul(delta, t2))
    c4 = add(add(mul(pk[0], u), e1), mul(delta, h))

    return e_enc, c0, c1, c2, c3, c4

def verifier(sk, proof):
    e_enc, c_0, c_1, c_2, c_3, c_4 = proof
    m_prime1 = np.poly1d(np.round(add(mul(e_enc, sk), c_0) * t / q) % t)

    print("m_prime1: ", m_prime1)
    m_prime2 = np.poly1d(np.round(add(mul(e_enc, sk), c_1) * t / q) % t)
    print("m_prime2: ", m_prime2)
    m_prime3 = np.poly1d(np.round(add(mul(e_enc, sk), c_2) * t / q) % t)
    print("m_prime3: ", m_prime3)
    m_prime4 = np.poly1d(np.round(add(mul(e_enc, sk), c_3) * t / q) % t)
    print("m_prime: ", m_prime4)
    m_prime5 = np.poly1d(np.round(add(mul(e_enc, sk), c_4) * t / q) % t)
    print("m_prime5: ", m_prime5)

    adec_alpha = mod_horner(m_prime1, alpha, t)
    bdec_alpha = mod_horner(m_prime2, alpha, t)
    cdec_alpha = mod_horner(m_prime3, alpha, t)
    tdec_alpha = mod_horner(m_prime4, alpha, t)
    hdec_alpha = mod_horner(m_prime5, alpha, t)
    print("adec_alpha: ", adec_alpha)
    print("bdec_alpha: ", bdec_alpha)
    print("cdec_alpha: ", cdec_alpha)
    print("tdec_alpha: ", tdec_alpha)
    print("hdec_alpha: ", hdec_alpha)

    left3 = (adec_alpha * bdec_alpha - cdec_alpha) % t
    right3 = (tdec_alpha * hdec_alpha) % t

    print("left3: ", left3)
    print("right3: ", right3)
    print("res3: ", left3 == right3)
    return (left3 == right3)

if __name__ == "__main__":
    alpha, sk, a2, e, pk, u, e1, e2 = setup()
    print(alpha)
    proof = prover(pk, u, e1, e2, alpha)
    verification = verifier(sk, proof)
    print(verification)
