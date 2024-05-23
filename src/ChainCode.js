'use strict';

const { Contract } = require('fabric-contract-api');

class SmartContract extends Contract {
    async VerifyProof(ctx, proof) {
        const [c_0, c_1] = proof;

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
    }

}

module.exports = SmartContract;
