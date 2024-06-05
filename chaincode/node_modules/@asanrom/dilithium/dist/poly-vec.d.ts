import { Polynomium } from "./poly";
export declare class PolynomiumVector {
    polynomiums: Polynomium[];
    static randomVec(rho: Uint8Array, eta: number, length: number, nonce: number): PolynomiumVector;
    static randomVecGamma1(seed: Uint8Array, length: number, gamma1: number, nonce: number): PolynomiumVector;
    static expandA(rho: Uint8Array, k: number, l: number): PolynomiumVector[];
    constructor(size: number);
    ntt(): PolynomiumVector;
    reduce(): void;
    decompose(gamma2: number): PolynomiumVector[];
    invnttTomont(): void;
    add(other: PolynomiumVector): PolynomiumVector;
    sub(other: PolynomiumVector): PolynomiumVector;
    caddq(): void;
    shift(): PolynomiumVector;
    powerRound(): PolynomiumVector[];
    pointwiseMontgomery(u: Polynomium): PolynomiumVector;
    private static pointwiseAccMontgomery;
    mulMatrixPointwiseMontgomery(M: PolynomiumVector[]): PolynomiumVector;
    size(): number;
    chknorm(bound: number): boolean;
    toString(): string;
}
