import galois
import numpy as np
import fromLWEtoR1CS as r1
import fromR1CStoQAP as qap
import sys
import timeit

order = 73
GF = galois.GF(order)
q = 8929
t = 73
d = 57
delta = q // t
p_q = np.poly1d([1] + ([0] * (d - 1)) + [1])

def polymod(poly, poly_mod, coeff_mod):
    return np.poly1d(np.floor(np.polydiv(poly, poly_mod)[1]) % coeff_mod)

def addition(poly_mod, coeff_mod):
    return lambda a, b: np.poly1d(polymod(np.polyadd(a, b), poly_mod, coeff_mod))

def multiplication(poly_mod, coeff_mod):
    return lambda a, b: np.poly1d(polymod(np.polymul(a, b), poly_mod, coeff_mod))

add = addition(p_q, q)
mul = multiplication(p_q, q)

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

    T = galois.Poly([1, order - 1], field=GF)
    for i in range(2, len(A) + 1):
        T = T * galois.Poly([1, order - i], field=GF)

    H = (Ua * Va - Wa) // T

    left = Ua(alpha) * Va(alpha) - Wa(alpha)
    right = H(alpha) * T(alpha)

    c1 = add(add(mul(pk[0], u), e1), mul(delta, left))
    c2 = add(add(mul(pk[0], u), e1), mul(delta, right))

    return c1, c2

def verifier(proof):
    c_0, c_1 = proof
    check = add(c_0, -c_1)
    return int(check.coefficients[0]) == 0

def main():
    s = np.array([1, 5, 8, 1, 8, 4, 7, 1, 0, 6, 4, 9, 1, 0, 8, 0, 1, 5, 0, 3, 2, 1, 0, 0, 1, 1, 1, 0, 1, 4, 0, 0, 0, 6, 0, 0, 1, 0, 0, 0, 1, 5, 0, 0, 2, 4, 4, 4, 6, 6, 7, 0, 0, 1, 5, 5, 7])
    alpha, sk, a2, e, pk, u, e1, e2 = setup()

    start_prover = timeit.default_timer()
    proof = prover(pk, u, e1, alpha, s)
    end_prover = timeit.default_timer()

    start_verifier = timeit.default_timer()
    verification = verifier(proof)
    end_verifier = timeit.default_timer()

    print("verification result:", verification)
    print("proof generation runtime: ", end_prover - start_prover, "s")
    print("verification time: ", end_verifier - start_verifier, "s")
    print("proof size: ", sys.getsizeof(proof), " bytes")

    return end_prover - start_prover, end_verifier - start_verifier

if __name__ == "__main__":
    main()
    prover_times = []
    verifier_times = []
    for _ in range(50):
        prover_time, verifier_time = main()
        prover_times.append(prover_time)
        verifier_times.append(verifier_time)

    average_prover_time = sum(prover_times) / len(prover_times)
    average_verifier_time = sum(verifier_times) / len(verifier_times)

    min_prover_time = min(prover_times)
    max_prover_time = max(prover_times)
    min_verifier_time = min(verifier_times)
    max_verifier_time = max(verifier_times)

    print(f"Average proof generation time over 50 runs: {average_prover_time} s")
    print(f"Lowest proof generation time over 50 runs: {min_prover_time} s")
    print(f"Highest proof generation time over 50 runs: {max_prover_time} s")

    print(f"Average verification time over 50 runs: {average_verifier_time} s")
    print(f"Lowest verification time over 50 runs: {min_verifier_time} s")
    print(f"Highest verification time over 50 runs: {max_verifier_time} s")
