import numpy as np
import sys
import lwe
import fromR1CStoQAP as qap
import galois

order = 2**8 + 1
GF = galois.GF(order)

np.set_printoptions(threshold=sys.maxsize)

n = 2
q = 257

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

def commit_witness(witness, public_key):
    # Containers for 'a' values and ciphertexts
    a_values = []
    ciphertexts = []

    # Iterate over each witness element
    for element in witness:
        # Encrypt each element using the LWE encryption function
        a, ciphertext = lwe.encrypt2(element, public_key)

        # Store the 'a' values and ciphertexts
        a_values.append(a)
        ciphertexts.append(int(round(ciphertext[0])))

    # Return the encryption components
    return a_values, ciphertexts


def homomorphic_add(enc_val1, enc_val2):
    print("1: ", enc_val1)
    print("2: ", enc_val2[1][0])
    # Homomorphic addition in LWE is straightforward: add the ciphertexts
    return int(enc_val1) + enc_val2[1][0]

def construct_encrypted_evaluations(As, Bs, Cs, encrypted_witness, alpha, public_key):
    # Using the homomorphic addition to simulate the evaluation of QAP polynomials
    enc_A_alpha = homomorphic_add(GF(encrypted_witness[0]), lwe.encrypt2(int(As(alpha)), public_key))
    enc_B_alpha = homomorphic_add(GF(encrypted_witness[1]), lwe.encrypt2(int(Bs(alpha)), public_key))
    enc_C_alpha = homomorphic_add(GF(encrypted_witness[2]), lwe.encrypt2(int(Cs(alpha)), public_key))
    return enc_A_alpha, enc_B_alpha, enc_C_alpha

def construct_proof(alpha, witness, T):
    # Evaluate A, B, C polynomials at alpha
    A, B, C, As, Bs, Cs = qap.polySum(witness)

    As_alpha = As(alpha)
    Bs_alpha = Bs(alpha)
    Cs_alpha = Cs(alpha)
    # Construct the solution polynomial H(x) from the witness and evaluate it at alpha
    H = (As * Bs - Cs) // T
    H_alpha = H(alpha)

    # QAP relation check at alpha
    left_side = np.sum(As_alpha * Bs_alpha) - np.sum(Cs_alpha)
    right_side = H_alpha * T(alpha)
    assert left_side == right_side, f"QAP relation does not hold: {left_side} != {right_side}"
    print("it passed!!")
    return {'H_alpha': H_alpha, 'A_alpha': As_alpha, 'B_alpha': Bs_alpha, 'C_alpha': Cs_alpha}



def main():
    #get the r1cs
    A, B, C = LWEToR1CS_transform()
    witness = np.array([1, 22, 38, 21, 14, 37, 19, 1, 0, 1, 1, 21, 0, 37, 0, 21, 37])
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
    a_values, ciphertexts = commit_witness(witness, public_key)
    print("witness: ", witness)
    print("encrypted witness: ", ciphertexts)
    alpha = np.random.randint(q)
    construct_proof(alpha, GF(witness), T)
    A, B, C, As, Bs, Cs = qap.polySum(GF(witness))
    construct_encrypted_evaluations(As, Bs, Cs, ciphertexts, alpha, public_key)
    # proof, witness = prover()




if __name__ == "__main__":
    main()