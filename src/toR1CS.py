import src.KeyGen as kg
import numpy as np

def generate_r1cs(a, t, u, v, m, s, e1, e2, q, m_to_check):
    n = len(a)
    r1cs = []

    # Constraint 1: t = a * s + e
    eq1 = [f"{t[i, 0]} - (" + " - ".join([f"{a[i, j]} * s_{j}" for j in range(n)]) + f" - e_{i})" for i in range(n)]
    r1cs.append(eq1)

    # Constraint 2: u = r * a + e1
    eq2 = [f"{u[i]} - (" + " - ".join([f"{a[i, j]} * r_{j}" for j in range(n)]) + f" - e1_{i})" for i in range(n)]
    r1cs.append(eq2)

    # Constraint 3: v = r * t + e2 + (m * ((q+1) >> 1))
    eq3 = [f"{v} - (" + " - ".join([f"{t[j, 0]} * r_{j}" for j in range(n)]) + f" - e2 - (m * {(q+1) >> 1}))"]
    r1cs.append(eq3)

    # Constraint 4: m_to_check = v - (u * s)
    eq4 = [f"{m_to_check} - (" + " - ".join([f"{u[j]} * s_{j}" for j in range(n)]) + f" - v)"]
    r1cs.append(eq4)

    # Constraint 5: Check correctness (m_to_check == m)
    eq5 = [f"{m_to_check} - {m[0]}"]
    r1cs.append(eq5)

    return r1cs

def main():
    n = 4
    q = 1666

    # Generate random parameters
    s = np.random.randint(3, size=(n, 1)) - 1
    a, t = kg.keyGen(s)
    m = np.random.randint(2, size=1)
    e1 = np.random.randint(3, size=n) - 1
    e2 = np.random.randint(3, size=1) - 1
    u, v = kg.encrypt(a, t, m, e1, e2)
    m_to_check = kg.decrypt(s, u, v)

    # Generate R1CS
    r1cs = generate_r1cs(a, t, u, v, m, s, e1, e2, q, m_to_check)

    # Display R1CS
    print("Generated R1CS:")
    for eq in r1cs:
        print(eq)

if __name__ == "__main__":
    main()
