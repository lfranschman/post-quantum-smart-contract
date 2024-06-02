import { DilithiumLevel, DilithiumParameterSpec } from "./param";
import { PolynomiumVector } from "./poly-vec";
import { DilithiumSignature } from "./signature";
/**
 * Private key parameters
 */
export interface DilithiumPrivateKeyParams {
    rho: Uint8Array;
    tr: Uint8Array;
    K: Uint8Array;
    s1: PolynomiumVector;
    s2: PolynomiumVector;
    t0: PolynomiumVector;
    s1Hat: PolynomiumVector;
    s2Hat: PolynomiumVector;
    t0Hat: PolynomiumVector;
    A: PolynomiumVector[];
}
/**
 * Private key
 */
export declare class DilithiumPrivateKey implements DilithiumPrivateKeyParams {
    rho: Uint8Array;
    tr: Uint8Array;
    K: Uint8Array;
    s1: PolynomiumVector;
    s2: PolynomiumVector;
    t0: PolynomiumVector;
    s1Hat: PolynomiumVector;
    s2Hat: PolynomiumVector;
    t0Hat: PolynomiumVector;
    A: PolynomiumVector[];
    spec: DilithiumParameterSpec;
    bytes: Uint8Array;
    /**
     * Parses a private key from a byte array
     * @param bytes The byte array
     * @param level The security level
     * @returns The private key
     */
    static fromBytes(bytes: Uint8Array, level: DilithiumLevel): DilithiumPrivateKey;
    /**
     * Parses private key from hex string
     * @param hex The string
     * @param level The security level
     * @returns The private key
     */
    static fromHex(hex: string, level: DilithiumLevel): DilithiumPrivateKey;
    /**
     * Parses private key from base 64 string
     * @param base64 The string
     * @param level The security level
     * @returns The private key
     */
    static fromBase64(base64: string, level: DilithiumLevel): DilithiumPrivateKey;
    constructor(spec: DilithiumParameterSpec, bytes: Uint8Array, params: DilithiumPrivateKeyParams);
    /**
     * Derives the public key from this private key
     * @returns The public key
     */
    derivePublicKey(): DilithiumPublicKey;
    /**
     * Creates a key pair from this key
     * @returns The key pair
     */
    toKeyPair(): DilithiumKeyPair;
    /**
     * @returns The private key as a byte array
     */
    getBytes(): Uint8Array;
    /**
     * @returns The private key as a hex string
     */
    toHex(): string;
    /**
     * @returns The private key as a base 64 string
     */
    toBase64(): string;
    /**
     * Signs a message
     * @param message The message
     * @returns The signature
     */
    sign(message: Uint8Array): DilithiumSignature;
}
/**
 * Public key
 */
export declare class DilithiumPublicKey {
    rho: Uint8Array;
    t1: PolynomiumVector;
    A: PolynomiumVector[];
    spec: DilithiumParameterSpec;
    bytes: Uint8Array;
    /**
     * Parses a public key from a byte array
     * @param bytes The byte array
     * @param level The level specification
     * @returns The public key
     */
    static fromBytes(bytes: Uint8Array, level: DilithiumLevel): DilithiumPublicKey;
    /**
     * Parses public key from hex string
     * @param hex The string
     * @param level The security level
     * @returns The public key
     */
    static fromHex(hex: string, level: DilithiumLevel): DilithiumPublicKey;
    /**
     * Parses public key from base 64 string
     * @param base64 The string
     * @param level The security level
     * @returns The public key
     */
    static fromBase64(base64: string, level: DilithiumLevel): DilithiumPublicKey;
    constructor(spec: DilithiumParameterSpec, bytes: Uint8Array, rho: Uint8Array, t1: PolynomiumVector, A: PolynomiumVector[]);
    /**
     * @returns The public key as a byte array
     */
    getBytes(): Uint8Array;
    /**
     * @returns The public key as a hex string
     */
    toHex(): string;
    /**
     * @returns The public key as a base 64 string
     */
    toBase64(): string;
    /**
     * Verifies a signature
     * @param message The message
     * @param sig The signature
     * @returns True only if the signature is valid
     */
    verifySignature(message: Uint8Array, sig: DilithiumSignature): boolean;
}
/**
 * Key pair (private + public)
 */
export declare class DilithiumKeyPair {
    private pub;
    private secret;
    /**
     * Creates a key pair from a private key
     * @param secret Private key
     * @returns The key pair
     */
    static fromPrivateKey(secret: DilithiumPrivateKey): DilithiumKeyPair;
    /**
     * Packs keys into a key pair
     * @param secret The private key
     * @param pub The public key
     * @returns
     */
    static fromKeys(secret: DilithiumPrivateKey, pub: DilithiumPublicKey): DilithiumKeyPair;
    /**
     * Generates a random keypair
     * @param level The security level
     * @param seed The seed for the generation of the key pair (if not prodiced, a random seed is generated)
     */
    static generate(level: DilithiumLevel, seed?: Uint8Array): DilithiumKeyPair;
    constructor(secret: DilithiumPrivateKey, pub?: DilithiumPublicKey);
    /**
     * @returns The private key
     */
    getPrivateKey(): DilithiumPrivateKey;
    /**
     * @returns The public key
     */
    getPublicKey(): DilithiumPublicKey;
    /**
     * Signs a message
     * @param message The message
     * @returns The signature
     */
    sign(message: Uint8Array): DilithiumSignature;
    /**
     * Verifies a signature
     * @param message The message
     * @param sig The signature
     * @returns True only if the signature is valid
     */
    verifySignature(message: Uint8Array, sig: DilithiumSignature): boolean;
}
