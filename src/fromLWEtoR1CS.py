import numpy as np
import sys

np.set_printoptions(threshold=sys.maxsize)

n = 2
q = 50

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

    # print(C_matrix)




    return A_matrix, B_matrix, C_matrix

def main():
    MA, MB, MC = LWEToR1CS_transform()
    s = np.array([1, 22, 38, 21, 14, 37, 19, 1, 0, 1, 1, 21, 0, 37, 0, 21, 37])

    print(np.matmul(MC, s) == (np.matmul(MA, s) * np.matmul(MB, s)))


if __name__ == "__main__":
    main()