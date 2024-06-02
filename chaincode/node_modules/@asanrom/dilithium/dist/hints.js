// Hints
"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.useHints = exports.makeHints = void 0;
const constants_1 = require("./constants");
const poly_1 = require("./poly");
const poly_vec_1 = require("./poly-vec");
function makeHints(gamma2, v0, v1) {
    const hints = {
        v: new poly_vec_1.PolynomiumVector(v0.size()),
        cnt: 0,
    };
    for (let i = 0; i < v0.size(); i++) {
        const hint = polyMakeHint(gamma2, v0.polynomiums[i], v1.polynomiums[i]);
        hints.cnt += hint.cnt;
        hints.v.polynomiums[i] = hint.v;
    }
    return hints;
}
exports.makeHints = makeHints;
function polyMakeHint(gamma2, a, b) {
    const hint = {
        v: new poly_1.Polynomium(constants_1.N),
        cnt: 0,
    };
    for (let i = 0; i < constants_1.N; i++) {
        hint.v.coef[i] = makeHint(gamma2, a.coef[i], b.coef[i]);
        hint.cnt += hint.v.coef[i];
    }
    return hint;
}
function makeHint(gamma2, a0, a1) {
    if (a0 <= gamma2 || a0 > constants_1.Q - gamma2 || (a0 === constants_1.Q - gamma2 && a1 === 0)) {
        return 0;
    }
    return 1;
}
function useHints(gamma2, u, h) {
    const res = new poly_vec_1.PolynomiumVector(u.size());
    for (let i = 0; i < res.size(); i++) {
        res.polynomiums[i] = polyUseHint(gamma2, u.polynomiums[i], h.polynomiums[i]);
    }
    return res;
}
exports.useHints = useHints;
function polyUseHint(gamma2, u, h) {
    const res = new poly_1.Polynomium(constants_1.N);
    for (let i = 0; i < constants_1.N; i++) {
        res.coef[i] = useHint(gamma2, u.coef[i], h.coef[i]);
    }
    return res;
}
function useHint(gamma2, a, hint) {
    let a0;
    let a1;
    a1 = (a + 127) >> 7;
    if (gamma2 === (constants_1.Q - 1) / 32) {
        a1 = (a1 * 1025 + (1 << 21)) >> 22;
        a1 &= 15;
    }
    else if (gamma2 === (constants_1.Q - 1) / 88) {
        a1 = (a1 * 11275 + (1 << 23)) >> 24;
        a1 ^= ((43 - a1) >> 31) & a1;
    }
    else {
        throw new Error("Invalid gamma2: " + gamma2);
    }
    a0 = a - a1 * 2 * gamma2;
    a0 -= (((constants_1.Q - 1) / 2 - a0) >> 31) & constants_1.Q;
    if (hint === 0) {
        return a1;
    }
    if (gamma2 === (constants_1.Q - 1) / 32) {
        if (a0 > 0)
            return (a1 + 1) & 15;
        else
            return (a1 - 1) & 15;
    }
    else if (gamma2 === (constants_1.Q - 1) / 88) {
        if (a0 > 0)
            return (a1 === 43) ? 0 : a1 + 1;
        else
            return (a1 === 0) ? 43 : a1 - 1;
    }
    else {
        throw new Error("Invalid gamma2: " + gamma2);
    }
}
