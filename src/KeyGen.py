import numpy as np
import galois

GF = galois.GF(257*257*257*257)
n = 4
q = 50


def mod_mult(m1, m2):
    return np.remainder(np.matmul(m1, m2), q)

def mod_add(m1, m2):
    return np.remainder((m1 + m2), q)

def mod_sub(m1, m2):
    return np.remainder((m1 - m2), q)

def keyGen(a, s, e):
    t = mod_mult(a, s)
    t = mod_add(t, e)
    return t

def encrypt(a, t, m, e1, e2):
    r = np.random.randint(3, size=n) - 1
    u = mod_mult(r, a)
    u = mod_add(u, e1)
    v = mod_mult(r, t)
    v = mod_add(v, e2)
    v = mod_add(v, (m*((q+1)>>1)))
    return u, v

def decrypt(s, u, v):
    f1 = mod_mult(u, s)
    m_with_error = mod_sub(v, f1)
    if (m_with_error[0] > (q//4)) and (m_with_error[0] < (3*(q//4))):
        return 1
    else:
        return 0

def main():
    a_first_column = np.random.randint(q, size=(n, 1))
    a = a_first_column
    for i in range(1, n):
        a = np.hstack((a, np.roll(a_first_column, i)))
    s = np.random.randint(3, size=(n, 1)) - 1
    e = np.random.randint(3, size=(n, 1)) - 1
    t = keyGen(a, s, e)
    print("Public Key (a):")
    print(a)
    print("\nother public Key (t):")
    print(t)
    m = np.random.randint(2, size=1)
    print("\nOriginal Message (m):")
    print(m)
    e1 = np.random.randint(3, size=n) - 1
    e2 = np.random.randint(3, size=1) - 1
    u, v = encrypt(a, t, m, e1, e2)
    m_to_check = decrypt(s, u, v)

    # Correct the comparison conditions


    print("\nDecrypted Message:")
    print(m_to_check)
    print("\nIs Decrypted Message equal to Original Message?")
    print(m_to_check == m[0])


def check():
    curve_order = 257
    L = GF(np.array([
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]]))
    R = GF(np.array([
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]))


    O = GF(np.array([
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 1, 0, 0, 0, 0, curve_order-1, curve_order-1, 0, 0, 0]]))

    # witness vector `s`
    s = GF(np.array([1, 169695, 3, 12, 7, 9, 27, 324, 144, 1728, 98]))

    # Calculate the result of the mutliplication of matrix, and check if it is equal to O*s
    result = np.matmul(O, s) == np.matmul(L, s) * np.matmul(R, s)

    def interpolate_column(col, nb):
        xs = GF(np.arange(1, nb + 1))
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

    ## computes all interpolated polynomials for L, R, and O
    U_polys = get_polys_of_matrix(L)
    V_polys = get_polys_of_matrix(R)
    W_polys = get_polys_of_matrix(O)

    def get_polys_of_S_matrix(matrix):
        polys = []
        nb_of_columns = len(matrix)
        # for each column
        for col_id in range(nb_of_columns):
            polys.append(galois.Poly([matrix[0]], field=GF))
        return np.array(polys)

    A_polys = get_polys_of_S_matrix(s)

    # Summing all the polynamials of the collection into one polynomial
    U = galois.Poly([0], field=GF)
    V = galois.Poly([0], field=GF)
    W = galois.Poly([0], field=GF)
    for i in range(len(U_polys)):
        U += U_polys[i]
        V += V_polys[i]
        W += W_polys[i]

    print("U: ", U)
    print("V: ", V)
    print("W: ", W)

    Ua = galois.Poly([0], field=GF)
    for i in range(len(s)):
        Ua += U_polys[i] * s[i]

    Va = galois.Poly([0], field=GF)
    for i in range(len(s)):
        Va += V_polys[i] * s[i]

    Wa = galois.Poly([0], field=GF)
    for i in range(len(s)):
        Wa += W_polys[i] * s[i]

    print("Ua: ", Ua)
    print("Va: ", Va)
    print("Wa: ", Wa)

if __name__ == "__main__":
    check()
