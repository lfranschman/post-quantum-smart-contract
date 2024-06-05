// Parameter specs
"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.DILITHIUM_LEVEL5_P = exports.DILITHIUM_LEVEL3_P = exports.DILITHIUM_LEVEL2_P = exports.DilithiumLevel = void 0;
const constants_1 = require("./constants");
const util_1 = require("./util");
/**
 * Dilithium security level
 */
class DilithiumLevel {
    /**
     * Get the definition for a security level of Dilithium
     * Level can be:
     *  - LEVEL 2
     *  - LEVEL 3
     *  - LEVEL 5
     * @param level The security level
     * @returns The specification and parameters for that level
     */
    static get(level) {
        let p;
        switch (level) {
            case 2:
                p = exports.DILITHIUM_LEVEL2_P;
                break;
            case 3:
                p = exports.DILITHIUM_LEVEL3_P;
                break;
            case 5:
                p = exports.DILITHIUM_LEVEL5_P;
                break;
            default:
                throw new Error("Invalid security level: " + level);
        }
        return new DilithiumLevel({
            level: level,
            rawParams: p,
            publicKeyLength: (0, util_1.getPublicKeyByteLength)(p),
            privateKeyLength: (0, util_1.getPrivateKeyByteLength)(p),
            signatureLength: (0, util_1.getSignatureByteLength)(p),
        });
    }
    constructor(spec) {
        this.spec = spec;
    }
}
exports.DilithiumLevel = DilithiumLevel;
exports.DILITHIUM_LEVEL2_P = {
    k: 4,
    l: 4,
    gamma1: 1 << 17,
    gamma2: Math.floor((constants_1.Q - 1) / 88),
    tau: 39,
    d: 13,
    chalentropy: 192,
    eta: 2,
    beta: 78,
    omega: 80,
};
exports.DILITHIUM_LEVEL3_P = {
    k: 6,
    l: 5,
    gamma1: 1 << 19,
    gamma2: Math.floor((constants_1.Q - 1) / 32),
    tau: 49,
    d: 13,
    chalentropy: 225,
    eta: 4,
    beta: 196,
    omega: 55,
};
exports.DILITHIUM_LEVEL5_P = {
    k: 8,
    l: 7,
    gamma1: 1 << 19,
    gamma2: Math.floor((constants_1.Q - 1) / 32),
    tau: 60,
    d: 13,
    chalentropy: 257,
    eta: 2,
    beta: 120,
    omega: 75,
};
