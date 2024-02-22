import numpy as np
import random
n = 4
k = 1
q = 13


def mod_mult(m1, m2):
    return (m1.dot(m2)) % q

def mod_add(m1, m2):
    return (m1 + m2) % q

def mod_sub(m1, m2):
    return (m1 - m2) % q

def keyGen(s):
    a = np.random.randint(q, size=(n, n))
    e = np.random.randint(3, size=(n, 1)) - 1
    t = mod_mult(a, s)
    t = mod_add(t, e)
    return a, t

def encrypt(a, t, m):
    r = np.random.randint(3, size=n) - 1
    e1 = np.random.randint(3, size=n) - 1
    e2 = np.random.randint(3, size=1) - 1
    u = mod_mult(r, a)
    u = mod_add(u, e1)
    v = mod_mult(r, t)
    v = mod_add(v, e2)
    v = mod_add(v, (m*((q+1)>>1)))
    return u, v

def decrypt(s, u, v):
    f1 = mod_mult(u, s)
    m_with_error = mod_sub(v, f1)
    return m_with_error

def main():
    s = np.random.randint(3, size=(n, 1)) - 1
    a, t = keyGen(s)
    print(a)
    print(t)
    m = np.random.randint(2, size=1)
    print(m)
    u, v = encrypt(a, t, m)
    me = decrypt(s, u, v)
    if ((me[0] > (q/4)) & (me < (3*(q/4)))):
        me[0] = 1
    else:
        me[0] = 0
    print(me[0])
    print(m[0])
    print(me[0] == m[0])

if __name__ == "__main__":
    main()