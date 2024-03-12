import galois
import numpy as np
import fromLWEtoR1CS as r1
import lwe as Kg
import fromR1CStoQAP as qap

order = 2**8 + 1
GF = galois.GF(order)
n = 2
q = 257

def get_crs():
    tau = np.random.randint(q)
    alpha = np.random.randint(q)
    beta = np.random.randint(q)
    gamma = np.random.randint(q)
    delta = np.random.randint(q)
    return tau, alpha, beta, gamma, delta


a_first_column = np.random.randint(q, size=(n, 1))
a = a_first_column
for i in range(1, n):
    a = np.hstack((a, np.roll(a_first_column, i)))
s = np.random.randint(3, size=(n, 1)) - 1
e = np.random.randint(3, size=(n, 1)) - 1

t = Kg.keyGen(a, s, e)

# n degree of polynomials
A, B, C = r1.LWEToR1CS_transform()

tau, alpha, beta, gamma, delta = get_crs()

n = len(A)-1
print("n check: ", n)
# public inputs (matrix a and vector t)
#l = 2
# private inputs = total inputs - public inputs
m = len(A[0])

T_coefficients = [1]
T = galois.Poly(T_coefficients, field=galois.GF(q))
for i in range(2, len(A) + 1):
    T *= galois.Poly([1, -i], field=galois.GF(q))

def getTau():
    # t(tau)
    tTau = T(tau)

    # powers of tau for [A] and [C] (on a)
    tau_a = []
    # powers of tau for [B] (on t)
    tau_t = []
    # Powers of tau for [A] and [C] (on lattice)
    for i in range(T.degree + 1):
        local_tau = int(tau ** i)
        tau_a.append(a * local_tau)
        tau_t.append(t * local_tau)

    # Powers of tau for h(tau)t(tau)
    t_a = []
    for i in range(T.degree + 1):
        local_tau = int(tau ** i)
        tmp = (local_tau * T(tau))
        t_a.append(a * int(tmp))
    return tau, tTau, tau_a, tau_t, t_a, a, t, s, e


tau, tTau, tau_A, tau_t, t_A, a, t, s, e = getTau()
# Print the results
# print("tau: ", tau)
# print("Ttau: ", tTau)
# print("Powers of tau for [A]:\n", tau_A)
# print("Powers of tau for [t]:\n", tau_t)
# print("Powers of tau for h(tau)t(tau) [A]:\n", t_A)
# print("Powers of tau for h(tau)t(tau) [t]:\n", t_t)


def hxBalancing(Ua, Va, Wa):
    H = (Ua * Va - Wa) // T
    print("t(x):", T)
    print("h(x):", H)
    # remainder should be null
    assert (Ua * Va - Wa) % T == 0, "Remainder not zero!"
    print("Ua(tau) * Va(tau) - Wa(tau) == H(T(tau)):", Ua(tau) * Va(tau) - Wa(tau) == H(tau) * T(tau))
    print("it is: ", Ua(tau) * Va(tau) - Wa(tau))
    print("deg(T) =", T.degree)
    print("deg(H) =", H.degree)
    return H, T


U, V, W, Ua, Va, Wa = qap.polySum()
hxBalancing(Ua, Va, Wa)


# import random
#
# # Parameters
# m = 10  # Number of constraints (adjust as needed)
#
# # Generate random values for Powers of Tau
# tau = [random.random() for _ in range(2 * m - 1)]
# alpha = random.random()
# beta = random.random()
#
# # Define generator points G1 and G2
# G1 = generate_g1()
# G2 = generate_g2()
#
# # Compute Powers of Tau
# tau_G1 = scalar_multiply(tau, G1)
# tau_G2 = scalar_multiply(tau[:m], G2)
# alpha_tau_G1 = scalar_multiply(tau[:m], G1, scalar=alpha)
# beta_tau_G1 = scalar_multiply(tau[:m], G1, scalar=beta)
# beta_G2 = scalar_multiply(G2, scalar=beta)
#
# # Store the Powers of Tau values for later use
#
# # Phase 2: Generate proving key
# # Assuming A, B, C are the summed polynomials as mentioned in the question
#
# gamma = random.random()
# delta = random.random()
#
# # Compute the polynomials L_i
# L = [(beta * A[i] + alpha * B[i] + C[i]) for i in range(m)]
#
# # Compute L_i(tau) * G1 using linear combinations
# L_tau_G1 = scalar_multiply(tau[:m], G1, scalar=delta**(-1), linear_combinations=L)
#
# # Compute proving key values
# proving_key_values = {
#     "alpha_beta_delta_G1": (alpha, beta, delta) * G1,
#     "tau_G1": tau_G1,
#     "delta_inv_L_tau_G1": delta**(-1) * L_tau_G1,
#     "tau_Zx_G1": scalar_multiply(tau[:m - 1], Zx_tau, G1),
#     "beta_delta_G2": (beta, delta) * G2,
#     "tau_G2": tau_G2
# }
#
# # Generate verification key
# verification_key_values = {
#     "alpha_G1": alpha * G1,
#     "gamma_inv_L_tau_G1": gamma**(-1) * L_tau_G1,
#     "beta_gamma_delta_G2": (beta, gamma, delta) * G2
# }

# Store the proving key and verification key values for later use


