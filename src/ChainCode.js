    'use strict';

    const { Contract, Context } = require('fabric-contract-api');

    class SmartContract extends Contract {

        async VerifyProof(ctx, proof) {
            // Implement logic to verify the proof
            const isValid = this.verifyProofLogic(proof);
            return isValid;
        }

        verifyProofLogic(proof) {
            [c_0, c_1] = proof
            check = add(c_0, -c_1)
            print("just a check: ", int(check.coefficients[0]))

            return int(check.coefficients[0]) == 0
        }
    }

    module.exports = SmartContract;
