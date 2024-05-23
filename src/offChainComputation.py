from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import numpy as np
import fhelweSNARK  # This module would contain the logic you've defined in your Python code
import galois
import logging

app = FastAPI()

class Witness(BaseModel):
    w1: list
    w2: list
    pk: list
    sk: list
    e: list
    u: list
    e1: list
    e2: list
    alpha: int
    t: int
    q: int

@app.post("/compute-proof")
async def compute_proof(witness: Witness):
    try:
        w1 = np.poly1d(witness.w1)
        w2 = np.poly1d(witness.w2)
        sk = np.poly1d(witness.sk)
        pk = [np.poly1d(witness.pk[0]), np.poly1d(witness.pk[1])]
        e = np.poly1d(witness.e)
        u = np.poly1d(witness.u)
        e1 = np.poly1d(witness.e1)
        e2 = np.poly1d(witness.e2)
        alpha = witness.alpha
        t = witness.t
        q = witness.q

        # Compute w
        w = np.array(np.round(fhelweSNARK.add(fhelweSNARK.mul(w2, sk), w1) * t / q) % t).astype(int)
        logging.info(f"w_check: {w}")

        proof = fhelweSNARK.prover(pk, u, e1, alpha, w)  # Function to compute the proof
        print("proof: ", proof)
        proof_data = [p.coeffs.tolist() for p in proof]

        return {"proof": proof_data}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000, loop="asyncio")
