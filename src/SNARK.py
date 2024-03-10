import numpy as np
import secrets as secrets
from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
import struct

SEED_BYTES = 16
q = 257 #GAMMA_P
n = 2 #GAMMA_N



class Crs:
    seed = secrets.token_bytes(SEED_BYTES)
    s = 0
    a_s = 0
    v = 0
    t = 0

class VerificationKey:
    alpha = np.random.randint(q)
    beta = np.random.randint(q)
    s = np.random.randint(q)
    sk = np.array([secrets.choice([-1, 0, 1]) for _ in range(n)]).reshape(-1, 1)




