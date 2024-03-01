import galois
import numpy as np
import fromLWEtoR1CS as r1


print("Initializing a large field...")
order = 2**8 + 1
GF = galois.GF(order)
print("Field initialized")
s = GF(np.array([1, 22, 38, 21, 14, 37, 19, 1, 0, 1, 1, 21, 0, 37, 0, 21, 37]))


# x = GF(np.array([1,2,3,4,5,6,7,8]))
# L3_column = GF(np.array([0,0,0,0,1,1,1,1]))
#
# L3_poly = galois.lagrange_poly(x, L3_column)
#
# print("The resulting polynomial is:\n", L3_poly, sep='')
#
# # Checking each column element
# print("L3_poly(1) == 0:", L3_poly(1) == 0)
# print("L3_poly(2) == 0:", L3_poly(2) == 0)
# print("L3_poly(3) == 0:", L3_poly(3) == 0)
# print("L3_poly(4) == 0:", L3_poly(4) == 0)
# print("L3_poly(5) == 1:", L3_poly(5) == 1)
# print("L3_poly(6) == 1:", L3_poly(6) == 1)
# print("L3_poly(7) == 1:", L3_poly(7) == 1)
# print("L3_poly(8) == 1:", L3_poly(8) == 1)


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




# print("Ua: ", Ua.degree)
# print("Va: ", Va.degree)
# print("Wa: ", Wa.degree)


