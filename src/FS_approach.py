import numpy as np

def Power2Round(r, D):
    # This function assumes r is a numpy array
    r0 = r % (2*D) // D
    return r - r0 * D, r0

def Decompose(r, gamma):
    r0 = r % gamma
    r1 = (r - r0) // gamma
    condition = (r1 == gamma - 1) & (r0 >= gamma // 2)
    r1[condition] = 0
    r0[condition] -= gamma
    return r1, r0

def HighBits(r, gamma):
    # This function assumes r is a numpy array
    return (r + gamma // 2) // gamma

def LowBits(r, gamma):
    return Decompose(r, gamma)[1]

def MakeGHint(z, r, gamma):
    # This function assumes z and r are numpy arrays
    m = (q - 1) // gamma
    r1 = HighBits(r, gamma)
    z1 = HighBits(z + r, gamma)
    return (z1 - r1) % (2*m)

def UseGHint(h, w1, A1z1, A2z2, c, tA1, gamma):
    # This function assumes h, w1, A1z1, A2z2, c, and tA1 are numpy arrays
    return (w1 - (h * A1z1 + c * A2z2 - c * tA1)) % (2*gamma)

# Given security parameters (for example purposes, they would be set appropriately for actual use)
q = 8380417  # Example modulus
gamma = (q - 1) // 16  # Example value for gamma

# Assuming D is defined somewhere in your code
D = q // 16  # Example D value for Power2Round, needs to be specified correctly
