import { Polynomium } from "./poly";
import { PolynomiumVector } from "./poly-vec";
export interface Hints {
    v: PolynomiumVector;
    cnt: number;
}
export interface Hint {
    v: Polynomium;
    cnt: number;
}
export declare function makeHints(gamma2: number, v0: PolynomiumVector, v1: PolynomiumVector): Hints;
export declare function useHints(gamma2: number, u: PolynomiumVector, h: PolynomiumVector): PolynomiumVector;
