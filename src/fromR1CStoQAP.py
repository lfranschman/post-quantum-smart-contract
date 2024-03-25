import galois
import numpy as np
import fromLWEtoR1CS as r1



print("Initializing a large field...")
order = 11
GF = galois.GF(order)
q = 2741
t = 11
d = 12
delta = q // t
# Polynomial modulus
p_q = np.poly1d([1] + ([0] * (d - 1)) + [1])
print("Field initialized")
# s = GF(np.array([1, 22, 38, 21, 14, 37, 19, 1, 0, 1, 1, 21, 0, 37, 0, 21, 37]))
s = GF(np.array([1, 2, 4, 1, 2, 3, 2, 1, 0, 1, 1, 1, 0, 3, 0, 1, 3]))
# s = GF(np.array([1, 6, 2, 1, 2, 3, 7, 0, 0, 1, 1, 2, 0, 3, 0, 1, 8]))
# x = GF(np.array([1,2,3,4,5,6,7,8]))


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

def interpolate_column(col, nb):
    xs = GF(np.arange(1, nb+1))
    # print(col)
    return galois.lagrange_poly(xs, col)

def get_polys_of_matrix(matrix):
    polys = []
    nb_of_rows = len(matrix)
    nb_of_columns = len(matrix[0])
    # for each column
    for col_id in range(nb_of_columns):
        column = []
        for row in range(nb_of_rows):
            column.append(matrix[row][col_id])
        polys.append(interpolate_column(GF(np.array(column)), nb_of_rows))
    return np.array(polys)

def polySum():
    A, B, C = r1.LWEToR1CS_transform()
    # witness = np.array([1, 2, 4, 1, 2, 3, 2, 1, 0, 1, 1, 1, 0, 3, 0, 1, 3])
    print(np.matmul(C, np.array(s)) == (np.matmul(A, np.array(s)) * np.matmul(B, np.array(s))))

    # print("Extra A check: ", A)
    ## computes all interpolated polynomials for L, R, and O
    U_polys = get_polys_of_matrix(A)
    V_polys = get_polys_of_matrix(B)
    W_polys = get_polys_of_matrix(C)

    # Summing all the polynamials of the collection into one polynomial
    U = galois.Poly([0], field=GF)
    V = galois.Poly([0], field=GF)
    W = galois.Poly([0], field=GF)
    # print("U first check: ", U_polys)
    # print("V first check: ", V_polys)
    # print("W first check: ", W_polys)
    for i in range(len(U_polys)):
        U += U_polys[i]
        V += V_polys[i]
        W += W_polys[i]
    # A_polys = get_polys_of_S_matrix(s)

    Ua = galois.Poly([0], field=GF)
    for i in range(len(s)):
        Ua += U_polys[i] * s[i]

    Va = galois.Poly([0], field=GF)
    for i in range(len(s)):
        Va += V_polys[i] * s[i]

    Wa = galois.Poly([0], field=GF)
    for i in range(len(s)):
        Wa += W_polys[i] * s[i]

    # print("U: ", U)
    # print("V: ", V)
    # print("W: ", W)
    #
    # print("Ua: ", Ua)
    # print("Va: ", Va)
    # print("Wa: ", Wa)
    #
    # print("U degree: ", U.degree)
    # print("V degree: ", V.degree)
    # print("W degree: ", W.degree)
    #
    # print("Ua degree: ", Ua.degree)
    # print("Va degree: ", Va.degree)
    # print("Wa degree: ", Wa.degree)

    return U, V, W, Ua, Va, Wa


def get_polys_of_S_matrix(matrix):
    polys = []
    nb_of_columns = len(matrix)
    # for each column
    for col_id in range(nb_of_columns):
        polys.append(galois.Poly([matrix[0]], field=GF))
    return np.array(polys)

def mod_horner(poly, x, modulus):
    result = 0
    for coefficient in poly.coefficients:
        result = (result * x + coefficient) % modulus
    return result

def to_poly1d(As, Bs, Cs, H, T):
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



if __name__ == "__main__":
    A, B, C = r1.LWEToR1CS_transform()

    U, V, W, Ua, Va, Wa = polySum()
    T = galois.Poly([1, order - 1], field=GF)
    ## Then multiply by the other values: (x-2)(x-3)...(x-7)
    for i in range(2, len(A) + 1):
        T = T * galois.Poly([1, order - i], field=GF)

    H = (Ua * Va - Wa) // T

    print("ok lets see: ", (Ua * Va - Wa) % T)
    print(Ua)
    print(Va)
    print(Wa)
    print(H)
    print(T)
    alpha = 6
    print(alpha)
    lhs = (Ua(alpha) * Va(alpha) - Wa(alpha))
    rhs = H(alpha) * T(alpha)
    print(lhs == rhs)

    a, b, c, t2, h = to_poly1d(Ua, Va, Wa, T, H)
    a_alpha = mod_horner(a, alpha, t)
    b_alpha = mod_horner(b, alpha, t)
    c_alpha = mod_horner(c, alpha, t)
    t_alpha = mod_horner(t2, alpha, t)
    h_alpha = mod_horner(h, alpha, t)
    print(a_alpha * b_alpha - c_alpha == t_alpha * h_alpha)

    add = addition(p_q, q)
    mul = multiplication(p_q, q)

    sk = np.poly1d(np.random.randint(0, 2, d))
    a2 = np.poly1d(np.random.randint(0, q, d) % q)
    e = np.poly1d(np.random.normal(0, 2, d).astype(int) % q)
    pk_0 = add(-mul(a2, sk), e)
    pk_1 = a2
    pk = (pk_0, pk_1)
    u = np.poly1d(np.random.randint(0, 2, d))
    e_1 = np.poly1d(np.random.normal(0, 2, d).astype(int) % q)
    e_2 = np.poly1d(np.random.normal(0, 2, d).astype(int) % q)

    e_enc = add(mul(pk[1], u), e_2)
    c_0 = add(add(mul(pk[0], u), e_1), mul(delta, a))
    c_1 = add(add(mul(pk[0], u), e_1), mul(delta, b))
    c_2 = add(add(mul(pk[0], u), e_1), mul(delta, c))
    c_3 = add(add(mul(pk[0], u), e_1), mul(delta, t2))
    c_4 = add(add(mul(pk[0], u), e_1), mul(delta, h))

    m_prime1 = np.poly1d(np.round(add(mul(e_enc, sk), c_0) * t / q) % t)
    print(a)
    print(m_prime1)
    m_prime2 = np.poly1d(np.round(add(mul(e_enc, sk), c_1) * t / q) % t)
    print(b)
    print(m_prime2)
    m_prime3 = np.poly1d(np.round(add(mul(e_enc, sk), c_2) * t / q) % t)
    print(c)
    print(m_prime3)
    m_prime4 = np.poly1d(np.round(add(mul(e_enc, sk), c_3) * t / q) % t)
    print(t2)
    print(m_prime4)
    m_prime5 = np.poly1d(np.round(add(mul(e_enc, sk), c_4) * t / q) % t)
    print(h)
    print(np.poly1d(np.round(add(mul(e_enc, sk), c_4) * t / q) % t))

    adec_alpha = mod_horner(m_prime1, alpha, t)
    bdec_alpha = mod_horner(m_prime2, alpha, t)
    cdec_alpha = mod_horner(m_prime3, alpha, t)
    tdec_alpha = mod_horner(m_prime4, alpha, t)
    hdec_alpha = mod_horner(m_prime5, alpha, t)
    print(adec_alpha * bdec_alpha - cdec_alpha == tdec_alpha * hdec_alpha)



