"""Boolean phase polynomials -- the algebra iq derives circuits from.

A phase polynomial is a multilinear polynomial over Boolean variables
(x^2 = x) whose value is a phase measured in *turns* (1 turn = 2*pi rad).
Because e^{2 pi i m} = 1 for every integer m, coefficients only matter
mod 1: any term whose coefficient is an integer *cancels*.  This is the
"cancelling parts" step that turns, e.g., the n^2 terms of jk/N in the
Fourier transform into the n(n+1)/2 gates of the QFT.
"""
from __future__ import annotations

import math
from fractions import Fraction
from itertools import combinations

import numpy as np

HALF = Fraction(1, 2)


def var_key(v: str):
    return (v[0], int(v[1:]))


def mono_str(m) -> str:
    return "·".join(sorted(m, key=var_key)) if m else "1"


def mod1(c):
    """Reduce to (-1/2, 1/2]."""
    r = c - math.floor(c)
    if r > HALF:
        r -= 1
    return r


def is_zero(c, tol=1e-9) -> bool:
    if isinstance(c, (int, Fraction)):
        return c == 0
    return abs(c) < tol


def nice(c, max_den=1 << 20, tol=1e-10):
    """Snap a float to a nearby rational (exact arithmetic thereafter)."""
    if isinstance(c, (int, Fraction)):
        return Fraction(c)
    f = Fraction(float(c)).limit_denominator(max_den)
    return f if abs(float(f) - c) < tol else float(c)


def fmt_coef(c) -> str:
    return str(c) if isinstance(c, Fraction) else f"{c:.6g}"


class Poly:
    """dict  frozenset(variables) -> coefficient (turns)."""

    def __init__(self, terms=None):
        self.terms = {}
        for m, c in (terms or {}).items():
            self._add(frozenset(m), c)

    def _add(self, m, c):
        v = self.terms.get(m, 0) + c
        if is_zero(v, 0):
            self.terms.pop(m, None)
        else:
            self.terms[m] = v

    # constructors
    @classmethod
    def var(cls, name):
        return cls({frozenset([name]): Fraction(1)})

    @classmethod
    def const(cls, c):
        return cls({frozenset(): Fraction(c)})

    @classmethod
    def register(cls, prefix: str, n: int):
        """Integer x = sum_a 2^a x_a as a polynomial in its bits."""
        return cls({frozenset([f"{prefix}{a}"]): Fraction(2 ** a) for a in range(n)})

    # arithmetic (Boolean: x*x = x)
    def __add__(self, other):
        other = other if isinstance(other, Poly) else Poly.const(other)
        p = Poly(self.terms)
        for m, c in other.terms.items():
            p._add(m, c)
        return p

    __radd__ = __add__

    def __neg__(self):
        return self * -1

    def __sub__(self, other):
        return self + (-other if isinstance(other, Poly) else -other)

    def __mul__(self, other):
        if not isinstance(other, Poly):
            return Poly({m: c * other for m, c in self.terms.items()})
        p = Poly()
        for m1, c1 in self.terms.items():
            for m2, c2 in other.terms.items():
                p._add(m1 | m2, c1 * c2)
        return p

    __rmul__ = __mul__

    def reduce_mod1(self, tol=1e-9):
        """Return (reduced poly, cancelled [(m, c)], kept [(m, c, c mod 1)])."""
        out, cancelled, kept = Poly(), [], []
        for m, c in self.terms.items():
            r = mod1(c)
            if is_zero(r, tol) or (not isinstance(r, Fraction) and abs(abs(r) - 1) < tol):
                cancelled.append((m, c))
            else:
                out.terms[m] = r
                kept.append((m, c, r))
        return out, cancelled, kept

    def sorted_terms(self):
        return sorted(self.terms.items(), key=lambda t: (len(t[0]), sorted(map(var_key, t[0]))))

    def degree(self):
        return max((len(m) for m in self.terms), default=0)

    def __eq__(self, other):
        return isinstance(other, Poly) and self.terms == other.terms

    def __str__(self):
        if not self.terms:
            return "0"
        return " + ".join(f"{fmt_coef(c)}·{mono_str(m)}" if m else fmt_coef(c)
                          for m, c in self.sorted_terms())


def moebius(table: np.ndarray, names) -> Poly:
    """Exact multilinear expansion of a Boolean function given as a (2,)*d table.

    f(x) = sum_S c_S prod_{i in S} x_i with c_S = sum_{T subset S} (-1)^{|S-T|} f(T).
    Integer coefficients => valid mod 1 for phases.
    """
    t = np.array(table, dtype=float).copy()
    d = t.ndim
    for ax in range(d):
        s1 = [slice(None)] * d
        s0 = [slice(None)] * d
        s1[ax], s0[ax] = 1, 0
        t[tuple(s1)] -= t[tuple(s0)]
    terms = {}
    for idx in np.ndindex(t.shape):
        c = t[idx]
        if abs(c) > 1e-12:
            terms[frozenset(names[i] for i, b in enumerate(idx) if b)] = nice(c)
    return Poly(terms)


def poly_of_transform(u: np.ndarray, n: int) -> Poly:
    """phi(j,k) with U[j,k] = e^{2 pi i phi(j,k)} / sqrt(N) (flat amplitudes)."""
    phi = np.angle(u * math.sqrt(2 ** n)) / (2 * math.pi)
    names = [f"j{n - 1 - i}" for i in range(n)] + [f"k{n - 1 - i}" for i in range(n)]
    return moebius(phi.reshape((2,) * (2 * n)), names)


def poly_of_diagonal(d: np.ndarray, n: int) -> Poly:
    """phi(k) with diag entry k = e^{2 pi i phi(k)}."""
    phi = np.angle(d) / (2 * math.pi)
    names = [f"k{n - 1 - i}" for i in range(n)]
    return moebius(phi.reshape((2,) * n), names)


def to_parity_basis(poly: Poly):
    """Rewrite monomials as sums of parities (XORs).

    prod_{i in S} x_i = 2^{1-|S|} sum_{T subset S, T != {}} (-1)^{|T|-1} XOR_{i in T} x_i
    Returns dict frozenset(T) -> coefficient (turns), constant dropped into key frozenset().
    """
    out = {}
    for m, c in poly.terms.items():
        if not m:
            out[m] = out.get(m, 0) + c
            continue
        s = sorted(m)
        scale = Fraction(1, 2 ** (len(s) - 1)) if isinstance(c, Fraction) else 2.0 ** (1 - len(s))
        for r in range(1, len(s) + 1):
            for t in combinations(s, r):
                out[frozenset(t)] = out.get(frozenset(t), 0) + c * scale * (-1) ** (r - 1)
    return out
