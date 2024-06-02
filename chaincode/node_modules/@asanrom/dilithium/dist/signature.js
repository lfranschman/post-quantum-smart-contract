// Signature algorithm
"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.DilithiumSignature = void 0;
const sha3_1 = require("sha3");
const constants_1 = require("./constants");
const hints_1 = require("./hints");
const poly_1 = require("./poly");
const poly_vec_1 = require("./poly-vec");
const util_1 = require("./util");
/**
 * Signature
 */
class DilithiumSignature {
    /**
     * Generates a signature
     * @param message The message to sign
     * @param privateKey The private key to use for the signature
     * @returns The signature
     */
    static generate(message, privateKey) {
        const spec = privateKey.spec;
        const signatureLength = (0, util_1.getSignatureByteLength)(spec);
        const sig = new Uint8Array(signatureLength);
        const A = privateKey.A;
        let conc = (0, util_1.mergeArrays)([privateKey.tr, message]);
        const mu = (0, util_1.crh)(conc);
        conc = (0, util_1.mergeArrays)([privateKey.K, mu]);
        const rhoprime = (0, util_1.crh)(conc);
        const s1 = privateKey.s1Hat;
        const s2 = privateKey.s2Hat;
        const t0 = privateKey.t0Hat;
        let kappa = 0;
        let finish = false;
        while (!finish) {
            const y = poly_vec_1.PolynomiumVector.randomVecGamma1(rhoprime, spec.l, spec.gamma1, kappa++);
            let z = y.ntt();
            const w = z.mulMatrixPointwiseMontgomery(A);
            w.reduce();
            w.invnttTomont();
            w.caddq();
            const res = w.decompose(spec.gamma2);
            (0, util_1.packw1)(spec.gamma2, res[1], sig);
            const s = new sha3_1.SHAKE(256);
            s.update(Buffer.from(mu));
            s.update(Buffer.from(sig.slice(0, res[1].size() * (0, util_1.getPolyW1PackedBytes)(spec.gamma2))));
            const bb = Buffer.alloc(constants_1.SEEDBYTES);
            s.digest({ buffer: bb, format: "binary" });
            sig.set(new Uint8Array(bb), 0);
            let cp = poly_1.Polynomium.generateChallenge(spec.tau, sig);
            cp = cp.ntt();
            z = s1.pointwiseMontgomery(cp);
            z.invnttTomont();
            z = z.add(y);
            z.reduce();
            if (z.chknorm(spec.gamma1 - spec.beta)) {
                continue;
            }
            let h = s2.pointwiseMontgomery(cp);
            h.invnttTomont();
            let w0 = res[0].sub(h);
            w0.reduce();
            if (w0.chknorm(spec.gamma2 - spec.beta)) {
                continue;
            }
            h = t0.pointwiseMontgomery(cp);
            h.invnttTomont();
            h.reduce();
            if (h.chknorm(spec.gamma2)) {
                continue;
            }
            w0 = w0.add(h);
            w0.caddq();
            const hints = (0, hints_1.makeHints)(spec.gamma2, w0, res[1]);
            if (hints.cnt > spec.omega) {
                continue;
            }
            (0, util_1.packSig)(spec.gamma1, spec.omega, sig, sig, z, hints.v);
            finish = true;
        }
        return new DilithiumSignature(spec, sig);
    }
    /**
     * Parses signature from byte array
     * @param bytes The byte array
     * @param level The level specification
     * @returns The signature
     */
    static fromBytes(bytes, level) {
        if (level.spec.signatureLength !== bytes.length) {
            throw new Error(`Invalid signature size. Expected ${level.spec.publicKeyLength} bytes, but found ${bytes.length} bytes`);
        }
        return new DilithiumSignature(level.spec.rawParams, bytes);
    }
    /**
     * Parses signature from hex string
     * @param hex The string
     * @param level The level specification
     * @returns The signature
     */
    static fromHex(hex, level) {
        return DilithiumSignature.fromBytes(new Uint8Array(Buffer.from(hex, "hex")), level);
    }
    /**
     * Parses signature from base 64 string
     * @param base64 The string
     * @param level The level specification
     * @returns The signature
     */
    static fromBase64(base64, level) {
        return DilithiumSignature.fromBytes(new Uint8Array(Buffer.from(base64, "base64")), level);
    }
    constructor(spec, bytes) {
        this.spec = spec;
        this.bytes = bytes;
    }
    /**
     * @returns The private key as a byte array
     */
    getBytes() {
        return this.bytes;
    }
    /**
     * @returns The signature as a hex string
     */
    toHex() {
        return Buffer.from(this.bytes).toString("hex");
    }
    /**
     * @returns The signature as a base 64 string
     */
    toBase64() {
        return Buffer.from(this.bytes).toString("base64");
    }
    /**
     * Verifies the signature against a message + public key
     * @param message The message
     * @param publicKey The public key
     * @returns True only if the signature is valid
     */
    verify(message, publicKey) {
        const spec = publicKey.spec;
        const signatureLength = (0, util_1.getSignatureByteLength)(spec);
        const sig = this.bytes;
        if (sig.length !== signatureLength) {
            return false; // Bad signature
        }
        let t1 = publicKey.t1;
        let off = 0;
        const c = sig.slice(0, constants_1.SEEDBYTES);
        off += constants_1.SEEDBYTES;
        let z = new poly_vec_1.PolynomiumVector(spec.l);
        for (let i = 0; i < spec.l; i++) {
            z.polynomiums[i] = poly_1.Polynomium.zunpack(spec.gamma1, sig, off);
            off += (0, util_1.getPolyZPackedBytes)(spec.gamma1);
        }
        const h = new poly_vec_1.PolynomiumVector(spec.k);
        let k = 0;
        for (let i = 0; i < h.size(); i++) {
            h.polynomiums[i] = new poly_1.Polynomium(constants_1.N);
            if ((sig[off + spec.omega + i] & 0xFF) < k || (sig[off + spec.omega + i] & 0xFF) > spec.omega) {
                return false;
            }
            for (let j = k; j < (sig[off + spec.omega + i] & 0xFF); j++) {
                /* Coefficients are ordered for strong unforgeability */
                if (j > k && (sig[off + j] & 0xFF) <= (sig[off + j - 1] & 0xFF)) {
                    return false;
                }
                h.polynomiums[i].coef[sig[off + j] & 0xFF] = 1;
            }
            k = (sig[off + spec.omega + i] & 0xFF);
        }
        for (let j = k; j < spec.omega; j++) {
            if (sig[off + j] !== 0) {
                return false;
            }
        }
        if (z.chknorm(spec.gamma1 - spec.beta)) {
            return false;
        }
        let mu = (0, util_1.crh)(publicKey.bytes);
        mu = (0, util_1.getSHAKE256Digest)(constants_1.CRHBYTES, mu, message);
        let cp = poly_1.Polynomium.generateChallenge(spec.tau, c);
        const A = publicKey.A;
        z = z.ntt();
        let w = z.mulMatrixPointwiseMontgomery(A);
        cp = cp.ntt();
        t1 = t1.shift();
        t1 = t1.ntt();
        t1 = t1.pointwiseMontgomery(cp);
        w = w.sub(t1);
        w.reduce();
        w.invnttTomont();
        w.caddq();
        w = (0, hints_1.useHints)(spec.gamma2, w, h);
        const buf = new Uint8Array((0, util_1.getPolyW1PackedBytes)(spec.gamma2) * w.size());
        (0, util_1.packw1)(spec.gamma2, w, buf);
        const c2 = (0, util_1.getSHAKE256Digest)(constants_1.SEEDBYTES, mu, buf);
        for (let i = 0; i < constants_1.SEEDBYTES; i++) {
            if (c[i] !== c2[i]) {
                return false;
            }
        }
        return true;
    }
}
exports.DilithiumSignature = DilithiumSignature;
