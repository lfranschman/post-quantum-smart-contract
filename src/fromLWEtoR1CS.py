import numpy as np

from src.KeyGen import keyGen

n = 4
q = 50

def LWEToR1CS_transform(a, s, e, t):

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

    # processing gates c0-cn^2
    for i in range(n*n):
        if (i % n) == 0:
            A_matrix[n*n + i][1 + 3*n + i + n*n] = 1
            A_matrix[n*n + i][1 + 3*n + n*n + i + 1] = 1

    i = 0
    while i < n:
        if ((i+1) % n) == 0:
              if (1 + 3 * n + i + 2 * n * n < num_variables):
                A_matrix[n*n + i][1 + 3 * n + i + 2 * n * n - 1] = 1
                A_matrix[n*n + i][1 + 2*n + i + n*n] = 1
                i = i + 1
    # for i in range(n * n):
    #     if ((i+1) % n) != 0 & (i % n) !=0:
    #         if(1 + 3 * n + i + 2* n * n < num_variables):
    #             A_matrix[n * n + i][1 + 3 * n + i + n * n + 1] = 1
    #             A_matrix[n * n + i][1 + 3 * n + i + 2 * n * n - 1] = 1






    # print(A_matrix[0])
    # print(A_matrix[1])
    # print(A_matrix[2])
    # print(A_matrix[3])
    # print(A_matrix[4])
    # print(A_matrix[5])
    # print(A_matrix[6])
    # print(A_matrix[7])
    # print(A_matrix[8])
    # print(A_matrix[9])
    # print(A_matrix[10])
    # print(A_matrix[11])
    # print(A_matrix[12])
    # print(A_matrix[13])
    # print(A_matrix[14])
    # print(A_matrix[15])
    print(A_matrix[16])
    print(A_matrix[17])
    print(A_matrix[18])
    print(A_matrix[19])
    print(A_matrix[20])
    print(A_matrix[21])
    print(A_matrix[22])
    print(A_matrix[23])
    print(A_matrix[24])
    print(A_matrix[25])
    print(A_matrix[26])
    print(A_matrix[27])
    print(A_matrix[28])
    print(A_matrix[29])
    print(A_matrix[30])
    print(A_matrix[31])


    # print(B_matrix[0])
    # print(B_matrix[1])
    # print(B_matrix[2])
    # print(B_matrix[3])
    # print(B_matrix[4])
    # print(B_matrix[5])
    # print(B_matrix[6])
    # print(B_matrix[7])
    # print(B_matrix[8])
    # print(B_matrix[9])
    # print(B_matrix[10])
    # print(B_matrix[11])
    # print(B_matrix[12])
    # print(B_matrix[13])
    # print(B_matrix[14])
    # print(B_matrix[15])
    # print(B_matrix[16])
    #
    # print(C_matrix[0])
    # print(C_matrix[1])
    # print(C_matrix[2])
    # print(C_matrix[3])
    # print(C_matrix[4])
    # print(C_matrix[5])
    # print(C_matrix[6])
    # print(C_matrix[7])
    # print(C_matrix[8])
    # print(C_matrix[9])
    # print(C_matrix[10])
    # print(C_matrix[11])
    # print(C_matrix[12])
    # print(C_matrix[13])
    # print(C_matrix[14])
    # print(C_matrix[15])
    # print(C_matrix[16])
    # print(C_matrix[17])

    return A_matrix

def main():
    a_first_column = np.random.randint(q, size=(n, 1))
    a = a_first_column
    for i in range(1, n):
        a = np.hstack((a, np.roll(a_first_column, i)))
    s = np.random.randint(3, size=(n, 1)) - 1
    t = keyGen(a, s)
    e = np.random.randint(3, size=(n, 1)) - 1

    LWEToR1CS_transform(a, s, e, t)
if __name__ == "__main__":
    main()