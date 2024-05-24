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

       
    proof = fhelweSNARK.prover(pk, u, e1, alpha, w)
    proof_response = [p.coeffs.tolist() for p in proof]

    
    print("proof_response: ", proof_response)
#    data = proof_response['proof']
    url = "http://localhost:2800/checkProof"
    payload = json.dumps(proof_response)
    headers = {
    'Content-Type': 'application/json'
    }
    response = requests.post(url, headers=headers, data=payload)
    logging.info(f"Check Proof response status: {response.status_code}")
    logging.info(f"Check Proof response content: {response.content}")
    

if __name__ == "__main__":
    asyncio.run(main())
