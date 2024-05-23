import requests
import logging
import json
import os
import numpy as np
import fhelweSNARK
import asyncio

# Configure basic logger
logging.basicConfig(level=logging.INFO)
order = 73
q = 8929
t = 73
d = 57
delta = q // t

# Polynomial modulus
p_q = np.poly1d([1] + ([0] * (d - 1)) + [1])

def send_witness_data(url, data):
    """
    Send witness data to the specified off-chain computation service and return the proof.
    """
    try:
        response = requests.post(url, json=data)
        response.raise_for_status()  # Will raise an HTTPError for bad responses (400 or 500 level responses)
        return response.json()
    except requests.exceptions.HTTPError as http_err:
        logging.error(f"HTTP error occurred: {http_err} - Status code: {response.status_code} - Response content: {response.content}")
    except Exception as err:
        logging.error(f"An error occurred: {err}")
    return None

async def main():
    off_chain_service_url = "http://localhost:8000/compute-proof"
    
    # Assuming fhelweSNARK.setup() returns these values correctly
    alpha, sk, a2, e, pk, u, e1, e2 = fhelweSNARK.setup()
    
    w = [1, 5, 8, 1, 8, 4, 7, 1, 0, 6, 4, 9, 1, 0, 8, 0, 1, 5, 0, 3, 2, 1, 0, 0, 1, 1, 1, 0, 1, 4, 0, 0, 0, 6, 0, 0, 1, 0, 0, 0, 1, 5, 0, 0, 2, 4, 4, 4, 6, 6, 7, 0, 0, 1, 5, 5, 7]

    w_encrypted_1 = fhelweSNARK.add(fhelweSNARK.add(fhelweSNARK.mul(pk[0], u), e1), fhelweSNARK.mul(delta, np.poly1d(w)))
    w_encrypted_2 = fhelweSNARK.add(fhelweSNARK.mul(pk[1], u), e2)

    logging.info(f"w_encrypted_1: {w_encrypted_1}")
    logging.info(f"w: {np.poly1d(w)}")
    logging.info(f"decryption_check: {np.array(np.round(fhelweSNARK.add(fhelweSNARK.mul(w_encrypted_2, sk), w_encrypted_1) * t / q) % t)}")
    
    witness_data = {
        "w1": w_encrypted_1.coeffs.tolist(),
        "w2": w_encrypted_2.coeffs.tolist(),
        "pk": (pk[0].coeffs.tolist(), pk[1].coeffs.tolist()),
        "sk": sk.coeffs.tolist(),
        "e": e.coeffs.tolist(),
        "u": u.coeffs.tolist(),
        "e1": e1.coeffs.tolist(),
        "e2": e2.coeffs.tolist(),
        "alpha": int(alpha),
        "t": t,
        "q": q        
    }

    proof_response = send_witness_data(off_chain_service_url, witness_data)

    if proof_response:
        logging.info("Success: Proof response received")
        data = proof_response['proof']
        url = "http://localhost:2800/checkProof"
        payload = json.dumps(data)
        headers = {
            'Content-Type': 'application/json'
        }
        response = requests.post(url, headers=headers, data=payload)
        logging.info(f"Check Proof response status: {response.status_code}")
        logging.info(f"Check Proof response content: {response.content}")
    else:
        logging.error("Failed to get valid response from service.")

if __name__ == "__main__":
    asyncio.run(main())
