// Utils
"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.packw1 = exports.packSig = exports.packPrvKey = exports.packPubKey = exports.getPrivateKeyByteLength = exports.getPublicKeyByteLength = exports.getSignatureByteLength = exports.getPolyZPackedBytes = exports.getPolyW1PackedBytes = exports.getPolyEtaPackedBytes = exports.crh = exports.getSHAKE256Digest = exports.clearBuffer = exports.mergeArrays = void 0;
const sha3_1 = require("sha3");
const constants_1 = require("./constants");
function mergeArrays(arrays) {
    let totalLength = 0;
    for (let i = 0; i < arrays.length; i++) {
        totalLength += arrays[i].length;
    }
    const merged = new Uint8Array(totalLength);
    let cp = 0;
    for (let i = 0; i < arrays.length; i++) {
        merged.set(arrays[i], cp);
        cp += arrays[i].length;
    }
    return merged;
}
exports.mergeArrays = mergeArrays;
function clearBuffer(buf) {
    for (let i = 0; i < buf.length; i++) {
        buf[i] = 0x00;
    }
}
exports.clearBuffer = clearBuffer;
function getSHAKE256Digest(sz, ...arr) {
    const data = mergeArrays(arr);
    const s = new sha3_1.SHAKE(256);
    s.update(Buffer.from(data));
    const out = Buffer.alloc(sz);
    s.digest({ buffer: out, format: "binary" });
    return new Uint8Array(out);
}
exports.getSHAKE256Digest = getSHAKE256Digest;
function crh(p) {
    return getSHAKE256Digest(constants_1.CRHBYTES, p);
}
exports.crh = crh;
function getPolyEtaPackedBytes(eta) {
    if (eta === 2) {
        return 96;
    }
    else if (eta === 4) {
        return 128;
    }
    else {
        throw new Error("Invalid etA: " + eta);
    }
}
exports.getPolyEtaPackedBytes = getPolyEtaPackedBytes;
function getPolyW1PackedBytes(gamma2) {
    let b;
    if (gamma2 === Math.floor((constants_1.Q - 1) / 88)) {
        b = 192;
    }
    else if (gamma2 === Math.floor((constants_1.Q - 1) / 32)) {
        b = 128;
    }
    else {
        throw new Error("Error invalid gamma2: " + gamma2);
    }
    return b;
}
exports.getPolyW1PackedBytes = getPolyW1PackedBytes;
function getPolyZPackedBytes(gamma1) {
    if (gamma1 === (1 << 17)) {
        return 576;
    }
    else if (gamma1 === (1 << 19)) {
        return 640;
    }
    else {
        throw new Error("Invalid gamma1: " + gamma1);
    }
}
exports.getPolyZPackedBytes = getPolyZPackedBytes;
function getSignatureByteLength(spec) {
    return (constants_1.SEEDBYTES +
        (spec.l * getPolyZPackedBytes(spec.gamma1)) +
        spec.omega +
        spec.k);
}
exports.getSignatureByteLength = getSignatureByteLength;
function getPublicKeyByteLength(spec) {
    return (constants_1.SEEDBYTES +
        (spec.k * constants_1.POLYT1_PACKEDBYTES));
}
exports.getPublicKeyByteLength = getPublicKeyByteLength;
function getPrivateKeyByteLength(spec) {
    let pkbytes;
    switch (spec.eta) {
        case 2:
            pkbytes = 96;
            break;
        case 4:
            pkbytes = 128;
            break;
        default:
            throw new Error("Illegal eta: " + spec.eta);
    }
    return ((2 * constants_1.SEEDBYTES) +
        constants_1.CRHBYTES +
        (spec.l * pkbytes) +
        (spec.k * pkbytes) +
        (spec.k * constants_1.POLYT0_PACKEDBYTES));
}
exports.getPrivateKeyByteLength = getPrivateKeyByteLength;
function packPubKey(rho, t) {
    const size = constants_1.SEEDBYTES + (t.size() * constants_1.POLYT1_PACKEDBYTES);
    const pk = new Uint8Array(size);
    for (let i = 0; i < constants_1.SEEDBYTES; i++) {
        pk[i] = rho[i];
    }
    const tl = t.size();
    for (let i = 0; i < tl; i++) {
        t.polynomiums[i].t1pack(pk, constants_1.SEEDBYTES + (i * constants_1.POLYT1_PACKEDBYTES));
    }
    return pk;
}
exports.packPubKey = packPubKey;
function packPrvKey(eta, rho, tr, K, t0, s1, s2) {
    let pkbytes;
    switch (eta) {
        case 2:
            pkbytes = 96;
            break;
        case 4:
            pkbytes = 128;
            break;
        default:
            throw new Error("Illegal eta: " + eta);
    }
    const keySize = ((2 * constants_1.SEEDBYTES) +
        constants_1.CRHBYTES +
        (s1.size() * pkbytes) +
        (s2.size() * pkbytes) +
        (s2.size() * constants_1.POLYT0_PACKEDBYTES));
    const buf = new Uint8Array(keySize);
    let off = 0;
    for (let i = 0; i < constants_1.SEEDBYTES; i++) {
        buf[off + i] = rho[i];
    }
    off += constants_1.SEEDBYTES;
    for (let i = 0; i < constants_1.SEEDBYTES; i++) {
        buf[off + i] = K[i];
    }
    off += constants_1.SEEDBYTES;
    for (let i = 0; i < constants_1.CRHBYTES; i++) {
        buf[off + i] = tr[i];
    }
    off += constants_1.CRHBYTES;
    for (let i = 0; i < s1.size(); i++) {
        s1.polynomiums[i].etapack(eta, buf, off);
        off += pkbytes;
    }
    for (let i = 0; i < s2.size(); i++) {
        s2.polynomiums[i].etapack(eta, buf, off);
        off += pkbytes;
    }
    for (let i = 0; i < t0.size(); i++) {
        t0.polynomiums[i].t0pack(buf, off);
        off += constants_1.POLYT0_PACKEDBYTES;
    }
    return buf;
}
exports.packPrvKey = packPrvKey;
function packSig(gamma1, omega, sig, c, z, h) {
    const pkBytes = getPolyZPackedBytes(gamma1);
    let off = 0;
    for (let i = 0; i < constants_1.SEEDBYTES; i++) {
        sig[i] = c[i];
    }
    off += constants_1.SEEDBYTES;
    const zlength = z.size();
    for (let i = 0; i < zlength; i++) {
        z.polynomiums[i].zpack(gamma1, sig, off);
        off += pkBytes;
    }
    /* Encode h */
    const hlength = h.size();
    for (let i = 0; i < omega + hlength; i++) {
        sig[off + i] = 0;
    }
    let k = 0;
    for (let i = 0; i < hlength; i++) {
        for (let j = 0; j < constants_1.N; j++) {
            if (h.polynomiums[i].coef[j] !== 0) {
                sig[off + k++] = j & 0xFF;
            }
        }
        sig[off + omega + i] = k & 0xFF;
    }
}
exports.packSig = packSig;
function packw1(gamma2, w, sig) {
    const pkBytes = getPolyW1PackedBytes(gamma2);
    let off = 0;
    const length = w.size();
    for (let i = 0; i < length; i++) {
        w.polynomiums[i].w1pack(gamma2, sig, off);
        off += pkBytes;
    }
}
exports.packw1 = packw1;
