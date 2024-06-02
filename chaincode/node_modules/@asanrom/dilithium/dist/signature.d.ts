import { DilithiumPrivateKey, DilithiumPublicKey } from "./key";
import { DilithiumLevel, DilithiumParameterSpec } from "./param";
/**
 * Signature
 */
export declare class DilithiumSignature {
    spec: DilithiumParameterSpec;
    bytes: Uint8Array;
    /**
     * Generates a signature
     * @param message The message to sign
     * @param privateKey The private key to use for the signature
     * @returns The signature
     */
    static generate(message: Uint8Array, privateKey: DilithiumPrivateKey): DilithiumSignature;
    /**
     * Parses signature from byte array
     * @param bytes The byte array
     * @param level The level specification
     * @returns The signature
     */
    static fromBytes(bytes: Uint8Array, level: DilithiumLevel): DilithiumSignature;
    /**
     * Parses signature from hex string
     * @param hex The string
     * @param level The level specification
     * @returns The signature
     */
    static fromHex(hex: string, level: DilithiumLevel): DilithiumSignature;
    /**
     * Parses signature from base 64 string
     * @param base64 The string
     * @param level The level specification
     * @returns The signature
     */
    static fromBase64(base64: string, level: DilithiumLevel): DilithiumSignature;
    constructor(spec: DilithiumParameterSpec, bytes: Uint8Array);
    /**
     * @returns The private key as a byte array
     */
    getBytes(): Uint8Array;
    /**
     * @returns The signature as a hex string
     */
    toHex(): string;
    /**
     * @returns The signature as a base 64 string
     */
    toBase64(): string;
    /**
     * Verifies the signature against a message + public key
     * @param message The message
     * @param publicKey The public key
     * @returns True only if the signature is valid
     */
    verify(message: Uint8Array, publicKey: DilithiumPublicKey): boolean;
}
