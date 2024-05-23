'use strict';

const { Contract } = require('fabric-contract-api');

class smartContract extends Contract {

    async VerifyProof(ctx, proofString) {
        try {
            // Parse the proofString into an array
            const proof = JSON.parse(proofString);
    
            // Ensure proof is an array with exactly two elements
            if (!Array.isArray(proof) || proof.length !== 2) {
                throw new Error('Invalid proof format');
            }
    
            const [c_0, c_1] = proof;
    
            // Ensure both elements are arrays
            if (!Array.isArray(c_0) || !Array.isArray(c_1)) {
                throw new Error('Invalid proof format');
            }
    
            // Check if the lengths of the arrays are equal
            if (c_0.length !== c_1.length) {
                return false;
            }
    
            // Verify each element in the arrays
            for (let i = 0; i < c_0.length; i++) {
                if (c_0[i] !== c_1[i]) {
                    return false;
                }
            }
    
            // If all elements match, return true
            return true;
        } catch (error) {
            console.error('Error in VerifyProof:', error);
            return false;
        }
    }


}

module.exports = smartContract;
