from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import numpy as np
import fhelweSNARK  # This module would contain the logic you've defined in your Python code
import galois

app = FastAPI()

class Witness(BaseModel):
    w: list  # This will be the structure of your input witness

@app.post("/compute-proof")
async def compute_proof(witness: Witness):
    try:
        w = np.array(witness.w)
        alpha, sk, a2, e, pk, u, e1, e2 = fhelweSNARK.setup()
        proof = fhelweSNARK.prover(pk, u, e1, alpha, w)  # Function to compute the proof
        print("proof: ", proof)
        proof_data = [p.coeffs.tolist() for p in proof]
        return {"proof": proof_data}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)