const { DilithiumKeyPair, DilithiumLevel, DilithiumPublicKey, DilithiumSignature } = require("@asanrom/dilithium");

function signProof(proof) {
    try {
        const level = DilithiumLevel.get(2); // Get the security level config (2, 3, or 5)
        const keyPair = DilithiumKeyPair.generate(level);
        const privateKey = keyPair.getPrivateKey();
        const publicKey = keyPair.getPublicKey();
        const message = new Uint8Array(Buffer.from(JSON.stringify(proof), "utf8")); // Removed extra parenthesis
        const signature = privateKey.sign(message).getBytes(); // Get raw byte data from signature
        return JSON.stringify({
            signature: Buffer.from(signature).toString('base64'),
            publicKey: Buffer.from(publicKey.getBytes()).toString('base64'),
            proof: proof // Including proof for convenience
        });
    } catch (error) {
        console.error("Error in signProof:", error);
        return JSON.stringify({ error: error.message });
    }
}

function verifySignature(publicKeyStr, proof, signatureStr) {
    try {
        const level = DilithiumLevel.get(2); // Ensure we use the same security level
        const publicKeyBytes = Uint8Array.from(Buffer.from(publicKeyStr, 'base64'));
        const publicKey = DilithiumPublicKey.fromBytes(publicKeyBytes, level);
        const message = new Uint8Array(Buffer.from(JSON.stringify(proof), "utf8"));
        const signatureBytes = Uint8Array.from(Buffer.from(signatureStr, 'base64'));
        const signature = DilithiumSignature.fromBytes(signatureBytes, level);
        return publicKey.verifySignature(message, signature);
    } catch (error) {
        console.error("Error in verifySignature:", error);
        return false;
    }
}

const [,, action, ...args] = process.argv;

if (action === "signProof") {
    const proof = JSON.parse(args[0]);
    const key_and_sig = signProof(proof);
    console.log(key_and_sig);
} else if (action === "verifySignature") {
    const [publicKeyStr, proofStr, signatureStr] = args;
    const proof = JSON.parse(proofStr);
    const valid = verifySignature(publicKeyStr, proof, signatureStr);
    console.log(valid);
}
