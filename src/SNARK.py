import numpy as np
import gmpy2
import lwe
from Crypto.PublicKey import RSA
from Crypto.Random import get_random_bytes
from Crypto.Cipher import PKCS1_OAEP

# Define constants used in the code
GAMMA_D = 256  # Example value, adjust according to actual needs
GAMMA_M = 256  # Example value, adjust according to actual needs
GAMMA_P = 2**16  # Example prime number, adjust as needed
CT_BYTES = 32  # Example value, adjust according to actual cryptographic details

class Proof:
    def __init__(self):
        self.h = None
        self.hat_h = None
        self.hat_v = None
        self.v_w = None
        self.b_w = None

    def clear(self):
        # Clearing is handled automatically by Python's garbage collector
        pass

class CRS:
    def __init__(self):
        self.seed = get_random_bytes(16)  # Adjust size as needed
        self.s = np.zeros((GAMMA_D, CT_BYTES), dtype=np.uint8)
        self.as_ = np.zeros((GAMMA_D, CT_BYTES), dtype=np.uint8)
        self.v = np.zeros((GAMMA_M, CT_BYTES), dtype=np.uint8)
        self.t = np.zeros(CT_BYTES, dtype=np.uint8)

    def clear(self):
        # Clearing is handled automatically by Python's garbage collector
        pass


def poly_evaluate(coefficients, x, modulus):
    """Evaluate polynomial at given point x with coefficients in a finite field defined by modulus."""
    result = 0
    for coeff in reversed(coefficients):
        result = (result * x + coeff) % modulus
    return result


def setup(crs, vrs, ssp):
    # Initialize the RNG with the seed from crs
    rng = np.random.default_rng(np.frombuffer(crs.seed, dtype=np.uint32))

    # # Define LWE parameters
    # lwe_n = GAMMA_D  # LWE dimension, adjust as needed
    # lwe_q = GAMMA_P  # LWE modulus, adjust as needed



    # Sample vrs attributes
    vrs.alpha = rng.integers(low=1, high=GAMMA_P)
    vrs.beta = rng.integers(low=1, high=GAMMA_P)
    vrs.s = rng.integers(low=1, high=GAMMA_P)

    # Generate LWE public key as a vector of random integers within the field defined by q
    vrs.sk = np.random.randint(low=0, high=GAMMA_P, size=GAMMA_D)

    # Encrypt using LWE and update crs attributes
    for i in range(GAMMA_D):
        # Encrypt s_i and as_i using LWE
        _, encrypted_s_i = lwe.encrypt2(s_i, vrs.sk)
        _, encrypted_as_i = lwe.encrypt2(as_i, vrs.sk)

        # Convert encrypted values to the required format and store in crs
        crs.s[i] = np.array(encrypted_s_i, dtype=np.uint8)[:CT_BYTES]
        crs.as_[i] = np.array(encrypted_as_i, dtype=np.uint8)[:CT_BYTES]

        # Update s_i and as_i for the next iteration
        s_i = (s_i * vrs.s) % GAMMA_P
        as_i = (as_i * vrs.s) % GAMMA_P

    # β * t(s)
    # Assuming ssp_t_offset points to the coefficients of the t(s) polynomial
    t_s_coeffs = ssp['t_s']  
    v_i_bs = (poly_evaluate(t_s_coeffs, vrs.s, GAMMA_P) * vrs.beta) % GAMMA_P
    _, encrypted_t_s = lwe.encrypt2(v_i_bs, vrs.sk)
    # Assuming crs.t can store the encrypted value directly
    crs.t = encrypted_t_s

    # β * v_i for each i
    for i in range(1, GAMMA_M):
        # Assuming ssp_v_offset function or mapping exists that gives the coefficients for v_i polynomials
        v_i_coeffs = ssp[f'v_{i}']  # Adjust based on your data structure
        v_i_bs = (poly_evaluate(v_i_coeffs, vrs.s, GAMMA_P) * vrs.beta) % GAMMA_P
        _, encrypted_v_i = lwe.encrypt2(v_i_bs, vrs.sk)
        # Store encrypted v_i in crs.v at the appropriate index
        crs.v[i - 1] = encrypted_v_i


