"""Target unitaries built directly from their mathematical definitions
(never from circuits), so derived circuits are checked against first principles."""
from __future__ import annotations

import math

import numpy as np


def qft(n, shift=0):
    """QFT|k> = N^{-1/2} sum_j e^{2 pi i j (k+shift) / N} |j>."""
    N = 2 ** n
    j = np.arange(N)[:, None]
    k = np.arange(N)[None, :]
    return np.exp(2j * math.pi * j * (k + shift) / N) / math.sqrt(N)


def iqft(n):
    return qft(n).conj().T


def hadamard(n):
    """H^{(x)n}: U[j,k] = (-1)^{j.k} / sqrt(N)."""
    N = 2 ** n
    j = np.arange(N)[:, None]
    k = np.arange(N)[None, :]
    par = np.vectorize(lambda x: bin(x).count("1") & 1)(j & k)
    return (1 - 2 * par) / math.sqrt(N)


def chirp(n):
    """Quadratic phase diag(e^{i pi k^2 / N})."""
    N = 2 ** n
    k = np.arange(N)
    return np.diag(np.exp(1j * math.pi * k * k / N))


def diffusion(n):
    """Grover diffusion 2|s><s| - I with |s> the uniform superposition."""
    N = 2 ** n
    s = np.full((N, 1), 1 / math.sqrt(N))
    return 2 * s @ s.T - np.eye(N)


def _perm(n, f):
    N = 2 ** n
    u = np.zeros((N, N))
    for x in range(N):
        u[f(x), x] = 1
    return u


def _bit(x, q):
    return (x >> q) & 1


def cx(n=2):  # control q0, target q1
    return _perm(2, lambda x: x ^ (_bit(x, 0) << 1))


def toffoli(n=3):  # controls q0,q1 target q2
    return _perm(3, lambda x: x ^ ((_bit(x, 0) & _bit(x, 1)) << 2))


def swap(n=2):
    return _perm(2, lambda x: (_bit(x, 0) << 1) | _bit(x, 1))


def fredkin(n=3):  # control q0 swaps q1,q2
    def f(x):
        if not _bit(x, 0):
            return x
        return (x & 1) | (_bit(x, 1) << 2) | (_bit(x, 2) << 1)
    return _perm(3, f)


def cz(n=2):
    return np.diag([1, 1, 1, -1]).astype(complex)


def ccz(n=3):
    d = np.ones(8, dtype=complex)
    d[7] = -1
    return np.diag(d)


def iswap(n=2):
    return np.array([[1, 0, 0, 0], [0, 0, 1j, 0], [0, 1j, 0, 0], [0, 0, 0, 1]])


def mcz(n):
    d = np.ones(2 ** n, dtype=complex)
    d[-1] = -1
    return np.diag(d)


# name -> (builder, default n, fixed size?)
TARGETS = {
    "qft": (qft, 3, False, "quantum Fourier transform"),
    "iqft": (iqft, 3, False, "inverse QFT"),
    "hadamard": (hadamard, 3, False, "Walsh-Hadamard transform H^n"),
    "chirp": (chirp, 3, False, "quadratic phase e^{iπk²/N}"),
    "diffusion": (diffusion, 3, False, "Grover diffusion 2|s><s|-I"),
    "mcz": (mcz, 3, False, "multi-controlled Z"),
    "cx": (cx, 2, True, "CNOT"),
    "cz": (cz, 2, True, "controlled-Z"),
    "swap": (swap, 2, True, "SWAP"),
    "iswap": (iswap, 2, True, "iSWAP"),
    "ccz": (ccz, 3, True, "doubly-controlled Z"),
    "toffoli": (toffoli, 3, True, "Toffoli / CCX"),
    "fredkin": (fredkin, 3, True, "Fredkin / controlled-SWAP"),
}


def build(name, n=None):
    f, dn, fixed, _ = TARGETS[name]
    n = dn if (n is None or fixed) else n
    return f(n), n
