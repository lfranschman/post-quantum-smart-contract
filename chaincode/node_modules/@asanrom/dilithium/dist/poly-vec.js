// Vector of Polynomiums
"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.PolynomiumVector = void 0;
const constants_1 = require("./constants");
const poly_1 = require("./poly");
class PolynomiumVector {
    static randomVec(rho, eta, length, nonce) {
        const pv = new PolynomiumVector(length);
        for (let i = 0; i < length; i++) {
            pv.polynomiums[i] = poly_1.Polynomium.genRandom(rho, eta, nonce++);
        }
        return pv;
    }
    static randomVecGamma1(seed, length, gamma1, nonce) {
        const pv = new PolynomiumVector(length);
        for (let i = 0; i < length; i++) {
            pv.polynomiums[i] = poly_1.Polynomium.genRandomGamma1(seed, length * nonce + i, constants_1.N, gamma1);
        }
        return pv;
    }
    static expandA(rho, k, l) {
        const A = [];
        for (let i = 0; i < k; i++) {
            const pv = new PolynomiumVector(l);
            for (let j = 0; j < l; j++) {
                pv.polynomiums[j] = poly_1.Polynomium.genUniformRandom(rho, (i << 8) + j);
            }
            A.push(pv);
        }
        return A;
    }
    constructor(size) {
        this.polynomiums = [];
        for (let i = 0; i < size; i++) {
            this.polynomiums.push(null);
        }
    }
    ntt() {
        const length = this.polynomiums.length;
        const pv = new PolynomiumVector(length);
        for (let i = 0; i < length; i++) {
            pv.polynomiums[i] = this.polynomiums[i].ntt();
        }
        return pv;
    }
    reduce() {
        const length = this.polynomiums.length;
        for (let i = 0; i < length; i++) {
            this.polynomiums[i].reduce();
        }
    }
    decompose(gamma2) {
        const length = this.polynomiums.length;
        const res = [
            new PolynomiumVector(length),
            new PolynomiumVector(length),
        ];
        for (let i = 0; i < length; i++) {
            const r = this.polynomiums[i].decompose(gamma2);
            res[0].polynomiums[i] = r[0];
            res[1].polynomiums[i] = r[1];
        }
        return res;
    }
    invnttTomont() {
        const length = this.polynomiums.length;
        for (let i = 0; i < length; i++) {
            this.polynomiums[i].invnttTomont();
        }
    }
    add(other) {
        const length = this.polynomiums.length;
        const pv = new PolynomiumVector(length);
        for (let i = 0; i < length; i++) {
            if (!other.polynomiums[i]) {
                continue;
            }
            pv.polynomiums[i] = this.polynomiums[i].add(other.polynomiums[i]);
        }
        return pv;
    }
    sub(other) {
        const length = this.polynomiums.length;
        const pv = new PolynomiumVector(length);
        for (let i = 0; i < length; i++) {
            if (!other.polynomiums[i]) {
                continue;
            }
            pv.polynomiums[i] = this.polynomiums[i].sub(other.polynomiums[i]);
        }
        return pv;
    }
    caddq() {
        const length = this.polynomiums.length;
        for (let i = 0; i < length; i++) {
            this.polynomiums[i].caddq();
        }
    }
    shift() {
        const length = this.polynomiums.length;
        const pv = new PolynomiumVector(length);
        for (let i = 0; i < length; i++) {
            pv.polynomiums[i] = this.polynomiums[i].shiftl();
        }
        return pv;
    }
    powerRound() {
        const length = this.polynomiums.length;
        const res = [
            new PolynomiumVector(length),
            new PolynomiumVector(length),
        ];
        for (let i = 0; i < length; i++) {
            const r = this.polynomiums[i].powerRound();
            res[0].polynomiums[i] = r[0];
            res[1].polynomiums[i] = r[1];
        }
        return res;
    }
    pointwiseMontgomery(u) {
        const length = this.polynomiums.length;
        const pv = new PolynomiumVector(length);
        for (let i = 0; i < length; i++) {
            pv.polynomiums[i] = u.pointwiseMontgomery(this.polynomiums[i]);
        }
        return pv;
    }
    static pointwiseAccMontgomery(u, v) {
        let w = u.polynomiums[0].pointwiseMontgomery(v.polynomiums[0]);
        const length = v.size();
        for (let i = 1; i < length; i++) {
            const t = u.polynomiums[i].pointwiseMontgomery(v.polynomiums[i]);
            w = w.add(t);
        }
        return w;
    }
    mulMatrixPointwiseMontgomery(M) {
        const pv = new PolynomiumVector(M.length);
        for (let i = 0; i < M.length; i++) {
            pv.polynomiums[i] = PolynomiumVector.pointwiseAccMontgomery(M[i], this);
        }
        return pv;
    }
    size() {
        return this.polynomiums.length;
    }
    chknorm(bound) {
        for (const poly of this.polynomiums) {
            if (poly.chknorm(bound)) {
                return true;
            }
        }
        return false;
    }
    toString() {
        return "{" + this.polynomiums.map(a => a.toString()).join(", ") + "}";
    }
}
exports.PolynomiumVector = PolynomiumVector;
