'use strict';

const { Contract } = require('fabric-contract-api');
const { DilithiumPublicKey, DilithiumSignature, DilithiumLevel } = require('@asanrom/dilithium');

class SmartContract extends Contract {

    async verifySignature(publicKeyStr, proof, signatureStr) {
        try {
            const level = DilithiumLevel.get(2); // Ensure we use the same security level
            const publicKeyBytes = Uint8Array.from(Buffer.from(publicKeyStr, 'base64'));
            const publicKey = DilithiumPublicKey.fromBytes(publicKeyBytes, level);
            const message = new Uint8Array(Buffer.from(JSON.stringify(proof), "utf8"));
            const signatureBytes = Uint8Array.from(Buffer.from(signatureStr, 'base64'));
            const signature = DilithiumSignature.fromBytes(signatureBytes, level);
            console.log('Message just before verify:', message);
            console.log('publicKey just before verify:', publicKey);
            console.log('signature just before verify:', signature);
		
            return publicKey.verifySignature(message, signature);
        } catch (error) {
            console.error("Error in verifySignature:", error);
            return false;
        }
    }




    async VerifyProof(ctx, proofString) {
        try {
            console.log('Received proofString:', proofString);
            // Parse the proofString into an object
            const payloadParsed = JSON.parse(proofString);
            const proof = payloadParsed.proof;
            const publicKey = payloadParsed.public_key; // Assume base64 encoding
            const signature = payloadParsed.signature; // Assume base64 encoding

            console.log('Parsed payload:', payloadParsed);
            console.log('received proof:', proof);
            console.log('received publicKey:', publicKey);
            console.log('received signature:', signature);

            // Ensure proof is correctly formatted
            if (!Array.isArray(proof)) {
                console.error('Invalid proof format: proof is not an array');
                throw new Error('Invalid proof format');
            }

            //Verify the signature
            const isSignatureValid = await this.verifySignature(publicKey, proof, signature);
            if (!isSignatureValid) {
                console.error('Signature verification failed');
                return false;
            }

            // If the signature is valid, return true
            return true;
        } catch (error) {
            console.error('Error in VerifyProof:', error);
            return false;
        }
    }
}

module.exports = SmartContract;
