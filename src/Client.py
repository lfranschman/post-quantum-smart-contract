import requests
import logging

# Configure basic logger
logging.basicConfig(level=logging.INFO)

def send_witness_data(url, data):
    """
    Send witness data to the specified off-chain computation service and return the proof.
    """
    try:
        response = requests.post(url, json=data)
        response.raise_for_status()  # Will raise an HTTPError for bad responses (400 or 500 level responses)
        return response.json()
    except requests.exceptions.HTTPError as http_err:
        logging.error(f"HTTP error occurred: {http_err} - Status code: {response.status_code}")
    except Exception as err:
        logging.error(f"An error occurred: {err}")
    return None

def main():
    # Endpoint of the off-chain service
    off_chain_service_url = "http://localhost:8000/compute-proof"

    # Your witness data in the format expected by your FastAPI app
    witness_data = {
        "w": [1, 2, 4, 1, 2, 3, 2, 1, 0, 1, 1, 1, 0, 3, 0, 1, 3]  # Example data
    }

    # Send witness data and get the proof
    proof_response = send_witness_data(off_chain_service_url, witness_data)

    if proof_response:
        print("Success:", proof_response)
        # Additional code to submit proof to blockchain could go here
    else:
        print("Failed to get valid response from service.")

if __name__ == "__main__":
    main()
