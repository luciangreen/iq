"""Gate definitions.

Conventions
-----------
* Basis index ``k = sum_b k_b * 2**b`` -- qubit ``b`` holds bit ``k_b``
  (little-endian, as in Qiskit).
* A k-qubit gate matrix is indexed with its *first listed* qubit as the most
  significant bit, so ``CX(c, t)`` has the textbook matrix with control first.
* Phase-type gates (P, CP, MCP and their named special cases Z, S, T, CZ, CCZ)
  multiply the amplitude by ``e^{i*theta}`` exactly when *all* their qubits are 1.
  They are symmetric in their qubits.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

import numpy as np

TAU = 2 * math.pi
_S2 = 1 / math.sqrt(2)


def _diag(*v):
    return np.diag(np.array(v, dtype=complex))


FIXED = {
    "I": np.eye(2, dtype=complex),
    "H": np.array([[1, 1], [1, -1]], dtype=complex) * _S2,
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": _diag(1, -1),
    "S": _diag(1, 1j),
    "Sdg": _diag(1, -1j),
    "T": _diag(1, np.exp(1j * math.pi / 4)),
    "Tdg": _diag(1, np.exp(-1j * math.pi / 4)),
    "CX": np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex),
    "CZ": _diag(1, 1, 1, -1),
    "SWAP": np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], dtype=complex),
    "CCZ": _diag(1, 1, 1, 1, 1, 1, 1, -1),
}
_ccx = np.eye(8, dtype=complex)
_ccx[[6, 7]] = _ccx[[7, 6]]
FIXED["CCX"] = _ccx

ARITY = {"I": 1, "H": 1, "X": 1, "Y": 1, "Z": 1, "S": 1, "Sdg": 1, "T": 1, "Tdg": 1,
         "P": 1, "RZ": 1, "CX": 2, "CZ": 2, "CP": 2, "SWAP": 2, "CCX": 3, "CCZ": 3}
PARAMETRIC = {"P", "RZ", "CP", "MCP"}
SYMMETRIC = {"CZ", "CP", "MCP", "SWAP", "CCZ"}
SELF_INVERSE = {"I", "H", "X", "Y", "Z", "CX", "CZ", "SWAP", "CCX", "CCZ"}
DAGGER = {"S": "Sdg", "Sdg": "S", "T": "Tdg", "Tdg": "T"}
# phase-type gates and their angle
NAMED_PHASE = {"Z": math.pi, "S": math.pi / 2, "Sdg": -math.pi / 2, "T": math.pi / 4,
               "Tdg": -math.pi / 4, "CZ": math.pi, "CCZ": math.pi}
DIAGONAL = set(NAMED_PHASE) | {"I", "P", "CP", "MCP", "RZ"}


@dataclass(frozen=True)
class Gate:
    name: str
    qubits: tuple
    param: float | None = None

    def __post_init__(self):
        q = tuple(int(x) for x in self.qubits)
        if len(set(q)) != len(q):
            raise ValueError(f"repeated qubit in {self.name}{q}")
        if self.name in SYMMETRIC:
            q = tuple(sorted(q))
        object.__setattr__(self, "qubits", q)
        if self.param is not None:
            object.__setattr__(self, "param", float(self.param))

    def __str__(self):
        p = f"({fmt_angle(self.param)})" if self.param is not None else ""
        return f"{self.name}{p}" + "".join(f" q{q}" for q in self.qubits)


def phase_gate(qubits, theta) -> Gate:
    """Gate multiplying by e^{i theta} when all ``qubits`` are 1."""
    qubits = tuple(qubits)
    name = {1: "P", 2: "CP"}.get(len(qubits), "MCP")
    return Gate(name, qubits, theta)


def phase_of(g: Gate):
    """(qubits, theta) if ``g`` is a phase-type gate, else None."""
    if g.name in NAMED_PHASE:
        return g.qubits, NAMED_PHASE[g.name]
    if g.name in ("P", "CP", "MCP"):
        return g.qubits, g.param
    return None


def gate_matrix(g: Gate) -> np.ndarray:
    if g.name in FIXED:
        return FIXED[g.name]
    if g.name in ("P", "CP", "MCP"):
        d = np.ones(2 ** len(g.qubits), dtype=complex)
        d[-1] = np.exp(1j * g.param)
        return np.diag(d)
    if g.name == "RZ":
        return _diag(np.exp(-0.5j * g.param), np.exp(0.5j * g.param))
    raise ValueError(f"unknown gate {g.name}")


def inverse(g: Gate) -> Gate:
    if g.name in SELF_INVERSE:
        return g
    if g.name in DAGGER:
        return Gate(DAGGER[g.name], g.qubits)
    if g.name in PARAMETRIC:
        return Gate(g.name, g.qubits, -g.param)
    raise ValueError(f"no inverse for {g.name}")


def is_diagonal(g: Gate) -> bool:
    return g.name in DIAGONAL


def norm_angle(theta: float) -> float:
    """Angle reduced to (-pi, pi]."""
    r = math.remainder(theta, TAU)
    if r <= -math.pi + 1e-12:
        r = math.pi
    return r


def fmt_angle(theta: float) -> str:
    """Pretty angle as a rational multiple of pi when possible."""
    f = Fraction(theta / math.pi).limit_denominator(1 << 16)
    if abs(float(f) * math.pi - theta) > 1e-9:
        return f"{theta:.6g}"
    if f == 0:
        return "0"
    num, den = f.numerator, f.denominator
    head = {1: "π", -1: "-π"}.get(num, f"{num}π")
    return head if den == 1 else f"{head}/{den}"


def parse_angle(s: str) -> float:
    """Parse '3pi/4', '-π/8', 'pi', '0.25' into radians."""
    s = s.strip().replace("π", "pi").replace("*", "")
    if "pi" not in s:
        return float(s)
    num, _, den = s.partition("/")
    coef = num.replace("pi", "")
    c = {"": 1.0, "+": 1.0, "-": -1.0}.get(coef)
    c = float(coef) if c is None else c
    return c * math.pi / (float(den) if den else 1.0)


def parse_gate_spec(spec: str):
    """'CP(pi/2)' -> ('CP', angle); 'H' -> ('H', None)."""
    spec = spec.strip()
    if "(" in spec:
        name, arg = spec[:-1].split("(", 1)
        return name, parse_angle(arg)
    return spec, None
