import numpy as np

n = 4
q = 50


def mod_mult(m1, m2):
    return np.dot(m1, m2) % q

def mod_add(m1, m2):
    return (m1 + m2) % q

def mod_sub(m1, m2):
    return (m1 - m2) % q

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

if __name__ == "__main__":
    main()
