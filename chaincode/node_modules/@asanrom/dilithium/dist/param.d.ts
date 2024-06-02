/**
 * Dilithium algorithm params
 */
export interface DilithiumParameterSpec {
    k: number;
    l: number;
    gamma1: number;
    gamma2: number;
    tau: number;
    d: number;
    chalentropy: number;
    eta: number;
    beta: number;
    omega: number;
}
/**
 * Dilithium security level (numeric)
 */
export type DilithiumLevelNumber = 2 | 3 | 5;
/**
 * Dilithium security level
 */
export declare class DilithiumLevel {
    /**
     * Get the definition for a security level of Dilithium
     * Level can be:
     *  - LEVEL 2
     *  - LEVEL 3
     *  - LEVEL 5
     * @param level The security level
     * @returns The specification and parameters for that level
     */
    static get(level: DilithiumLevelNumber): DilithiumLevel;
    spec: DilithiumLevelSpec;
    constructor(spec: DilithiumLevelSpec);
}
/**
 * Dilithium level specification
 */
export interface DilithiumLevelSpec {
    /**
     * Level
     */
    level: DilithiumLevelNumber;
    /**
     * Raw algorithm parameters
     */
    rawParams: DilithiumParameterSpec;
    /**
     * Public key length in bytes
     */
    publicKeyLength: number;
    /**
     * Private key length in bytes
     */
    privateKeyLength: number;
    /**
     * Signature length in bytes
     */
    signatureLength: number;
}
export declare const DILITHIUM_LEVEL2_P: DilithiumParameterSpec;
export declare const DILITHIUM_LEVEL3_P: DilithiumParameterSpec;
export declare const DILITHIUM_LEVEL5_P: DilithiumParameterSpec;
