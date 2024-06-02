// Dilithium key-pair
"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.DilithiumKeyPair = exports.DilithiumPublicKey = exports.DilithiumPrivateKey = void 0;
const constants_1 = require("./constants");
const randombytes_1 = __importDefault(require("randombytes"));
const poly_1 = require("./poly");
const poly_vec_1 = require("./poly-vec");
const util_1 = require("./util");
const signature_1 = require("./signature");
/**
 * Private key
 */
class DilithiumPrivateKey {
    /**
     * Parses a private key from a byte array
     * @param bytes The byte array
     * @param level The security level
     * @returns The private key
     */
    static fromBytes(bytes, level) {
        if (level.spec.privateKeyLength !== bytes.length) {
            throw new Error(`Invalid private key size. Expected ${level.spec.publicKeyLength} bytes, but found ${bytes.length} bytes`);
        }
        const parameterSpec = level.spec.rawParams;
        const POLYETA_PACKEDBYTES = (0, util_1.getPolyEtaPackedBytes)(parameterSpec.eta);
        let off = 0;
        const rho = new Uint8Array(constants_1.SEEDBYTES);
        for (let i = 0; i < constants_1.SEEDBYTES; i++) {
            rho[i] = bytes[i];
        }
        off += constants_1.SEEDBYTES;
        const key = new Uint8Array(constants_1.SEEDBYTES);
        for (let i = 0; i < constants_1.SEEDBYTES; i++) {
            key[i] = bytes[off + i];
        }
        off += constants_1.SEEDBYTES;
        const tr = new Uint8Array(constants_1.CRHBYTES);
        for (let i = 0; i < constants_1.CRHBYTES; i++) {
            tr[i] = bytes[off + i];
        }
        off += constants_1.CRHBYTES;
        const s1 = new poly_vec_1.PolynomiumVector(parameterSpec.l);
        for (let i = 0; i < parameterSpec.l; i++) {
            s1.polynomiums[i] = poly_1.Polynomium.etaunpack(parameterSpec.eta, bytes, off);
            off += POLYETA_PACKEDBYTES;
        }
        const s2 = new poly_vec_1.PolynomiumVector(parameterSpec.k);
        for (let i = 0; i < parameterSpec.k; i++) {
            s2.polynomiums[i] = poly_1.Polynomium.etaunpack(parameterSpec.eta, bytes, off);
            off += POLYETA_PACKEDBYTES;
        }
        const t0 = new poly_vec_1.PolynomiumVector(parameterSpec.k);
        for (let i = 0; i < parameterSpec.k; i++) {
            t0.polynomiums[i] = poly_1.Polynomium.t0unpack(bytes, off);
            off += constants_1.POLYT0_PACKEDBYTES;
        }
        // Precompute A, s0, s1 & t0hat
        const A = poly_vec_1.PolynomiumVector.expandA(rho, parameterSpec.k, parameterSpec.l);
        const s1Hat = s1.ntt();
        const s2Hat = s2.ntt();
        const t0Hat = t0.ntt();
        return new DilithiumPrivateKey(parameterSpec, bytes, {
            rho,
            K: key,
            tr,
            s1,
            s2,
            t0,
            A,
            s1Hat,
            s2Hat,
            t0Hat,
        });
    }
    /**
     * Parses private key from hex string
     * @param hex The string
     * @param level The security level
     * @returns The private key
     */
    static fromHex(hex, level) {
        return DilithiumPrivateKey.fromBytes(new Uint8Array(Buffer.from(hex, "hex")), level);
    }
    /**
     * Parses private key from base 64 string
     * @param base64 The string
     * @param level The security level
     * @returns The private key
     */
    static fromBase64(base64, level) {
        return DilithiumPrivateKey.fromBytes(new Uint8Array(Buffer.from(base64, "base64")), level);
    }
    constructor(spec, bytes, params) {
        this.spec = spec;
        this.bytes = bytes;
        this.rho = params.rho;
        this.tr = params.tr;
        this.K = params.K;
        this.s1 = params.s1;
        this.s2 = params.s2;
        this.t0 = params.t0;
        this.s1Hat = params.s1Hat;
        this.s2Hat = params.s2Hat;
        this.t0Hat = params.t0Hat;
        this.A = params.A;
    }
    /**
     * Derives the public key from this private key
     * @returns The public key
     */
    derivePublicKey() {
        const s1hat = this.s1.ntt();
        let t1 = s1hat.mulMatrixPointwiseMontgomery(this.A);
        t1.reduce();
        t1.invnttTomont();
        t1 = t1.add(this.s2);
        t1.caddq();
        const res = t1.powerRound();
        const pubbytes = (0, util_1.packPubKey)(this.rho, res[1]);
        return new DilithiumPublicKey(this.spec, pubbytes, this.rho, res[1], this.A);
    }
    /**
     * Creates a key pair from this key
     * @returns The key pair
     */
    toKeyPair() {
        return new DilithiumKeyPair(this);
    }
    /**
     * @returns The private key as a byte array
     */
    getBytes() {
        return this.bytes;
    }
    /**
     * @returns The private key as a hex string
     */
    toHex() {
        return Buffer.from(this.bytes).toString("hex");
    }
    /**
     * @returns The private key as a base 64 string
     */
    toBase64() {
        return Buffer.from(this.bytes).toString("base64");
    }
    /**
     * Signs a message
     * @param message The message
     * @returns The signature
     */
    sign(message) {
        return signature_1.DilithiumSignature.generate(message, this);
    }
}
exports.DilithiumPrivateKey = DilithiumPrivateKey;
/**
 * Public key
 */
class DilithiumPublicKey {
    /**
     * Parses a public key from a byte array
     * @param bytes The byte array
     * @param level The level specification
     * @returns The public key
     */
    static fromBytes(bytes, level) {
        if (level.spec.publicKeyLength !== bytes.length) {
            throw new Error(`Invalid public key size. Expected ${level.spec.publicKeyLength} bytes, but found ${bytes.length} bytes`);
        }
        const parameterSpec = level.spec.rawParams;
        let off = 0;
        const rho = new Uint8Array(constants_1.SEEDBYTES);
        for (let i = 0; i < constants_1.SEEDBYTES; i++) {
            rho[i] = bytes[i];
        }
        off += constants_1.SEEDBYTES;
        const p = new poly_vec_1.PolynomiumVector(parameterSpec.k);
        for (let i = 0; i < parameterSpec.k; i++) {
            p.polynomiums[i] = poly_1.Polynomium.t1unpack(bytes, off);
            off += constants_1.POLYT1_PACKEDBYTES;
        }
        // Precompute A
        const A = poly_vec_1.PolynomiumVector.expandA(rho, parameterSpec.k, parameterSpec.l);
        return new DilithiumPublicKey(parameterSpec, bytes, rho, p, A);
    }
    /**
     * Parses public key from hex string
     * @param hex The string
     * @param level The security level
     * @returns The public key
     */
    static fromHex(hex, level) {
        return DilithiumPublicKey.fromBytes(new Uint8Array(Buffer.from(hex, "hex")), level);
    }
    /**
     * Parses public key from base 64 string
     * @param base64 The string
     * @param level The security level
     * @returns The public key
     */
    static fromBase64(base64, level) {
        return DilithiumPublicKey.fromBytes(new Uint8Array(Buffer.from(base64, "base64")), level);
    }
    constructor(spec, bytes, rho, t1, A) {
        this.spec = spec;
        this.bytes = bytes;
        this.rho = rho;
        this.t1 = t1;
        this.A = A;
    }
    /**
     * @returns The public key as a byte array
     */
    getBytes() {
        return this.bytes;
    }
    /**
     * @returns The public key as a hex string
     */
    toHex() {
        return Buffer.from(this.bytes).toString("hex");
    }
    /**
     * @returns The public key as a base 64 string
     */
    toBase64() {
        return Buffer.from(this.bytes).toString("base64");
    }
    /**
     * Verifies a signature
     * @param message The message
     * @param sig The signature
     * @returns True only if the signature is valid
     */
    verifySignature(message, sig) {
        return sig.verify(message, this);
    }
}
exports.DilithiumPublicKey = DilithiumPublicKey;
/**
 * Key pair (private + public)
 */
class DilithiumKeyPair {
    /**
     * Creates a key pair from a private key
     * @param secret Private key
     * @returns The key pair
     */
    static fromPrivateKey(secret) {
        return new DilithiumKeyPair(secret);
    }
    /**
     * Packs keys into a key pair
     * @param secret The private key
     * @param pub The public key
     * @returns
     */
    static fromKeys(secret, pub) {
        return new DilithiumKeyPair(secret, pub);
    }
    /**
     * Generates a random keypair
     * @param level The security level
     * @param seed The seed for the generation of the key pair (if not prodiced, a random seed is generated)
     */
    static generate(level, seed) {
        const spec = level.spec.rawParams;
        let zeta;
        if (seed) {
            zeta = seed;
        }
        else {
            zeta = new Uint8Array((0, randombytes_1.default)(32));
        }
        const o = (0, util_1.getSHAKE256Digest)(3 * 32, zeta);
        const rho = o.slice(0, 32);
        const sigma = o.slice(32, 64);
        const K = o.slice(64, 96);
        const s1 = poly_vec_1.PolynomiumVector.randomVec(sigma, spec.eta, spec.l, 0);
        const s2 = poly_vec_1.PolynomiumVector.randomVec(sigma, spec.eta, spec.k, spec.l);
        const A = poly_vec_1.PolynomiumVector.expandA(rho, spec.k, spec.l);
        const s1Hat = s1.ntt();
        let t1 = s1Hat.mulMatrixPointwiseMontgomery(A);
        t1.reduce();
        t1.invnttTomont();
        t1 = t1.add(s2);
        t1.caddq();
        const res = t1.powerRound();
        const pubbytes = (0, util_1.packPubKey)(rho, res[1]);
        const tr = (0, util_1.crh)(pubbytes);
        const prvbytes = (0, util_1.packPrvKey)(spec.eta, rho, tr, K, res[0], s1, s2);
        const s2Hat = s2.ntt();
        const t0Hat = res[0].ntt();
        const privateKey = new DilithiumPrivateKey(spec, prvbytes, {
            rho,
            K,
            tr,
            s1,
            s2,
            t0: res[0],
            A,
            s1Hat,
            s2Hat,
            t0Hat,
        });
        const publicKey = new DilithiumPublicKey(spec, pubbytes, rho, res[1], A);
        return new DilithiumKeyPair(privateKey, publicKey);
    }
    constructor(secret, pub) {
        this.secret = secret;
        this.pub = pub || secret.derivePublicKey();
    }
    /**
     * @returns The private key
     */
    getPrivateKey() {
        return this.secret;
    }
    /**
     * @returns The public key
     */
    getPublicKey() {
        return this.pub;
    }
    /**
     * Signs a message
     * @param message The message
     * @returns The signature
     */
    sign(message) {
        return this.secret.sign(message);
    }
    /**
     * Verifies a signature
     * @param message The message
     * @param sig The signature
     * @returns True only if the signature is valid
     */
    verifySignature(message, sig) {
        return this.pub.verifySignature(message, sig);
    }
}
exports.DilithiumKeyPair = DilithiumKeyPair;
