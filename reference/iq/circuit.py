"""Quantum circuits: unitary simulation, equivalence, OpenQASM 3 and ASCII drawing."""
from __future__ import annotations

import math
from collections import Counter

import numpy as np

from .gates import Gate, gate_matrix, fmt_angle, phase_of

QASM_NAMES = {"H": "h", "X": "x", "Y": "y", "Z": "z", "S": "s", "Sdg": "sdg", "T": "t",
              "Tdg": "tdg", "P": "p", "RZ": "rz", "CX": "cx", "CZ": "cz", "CP": "cp",
              "SWAP": "swap", "CCX": "ccx", "I": "id"}


def apply_gate(tensor: np.ndarray, g: Gate, n: int) -> np.ndarray:
    """Apply gate to a tensor of shape (2,)*n + rest (axis n-1-q is qubit q)."""
    axes = [n - 1 - q for q in g.qubits]
    k = len(axes)
    m = gate_matrix(g)
    t = np.moveaxis(tensor, axes, list(range(k)))
    shape = t.shape
    t = (m @ t.reshape(2 ** k, -1)).reshape(shape)
    return np.moveaxis(t, list(range(k)), axes)


class Circuit:
    def __init__(self, n: int, gates=(), global_phase: float = 0.0, output_map=None):
        self.n = n
        self.gates: list[Gate] = list(gates)
        self.global_phase = global_phase
        # output_map[a] = qubit holding output bit a (None means identity mapping)
        self.output_map = output_map

    # -- building -----------------------------------------------------------
    def add(self, name, qubits, param=None):
        if isinstance(qubits, int):
            qubits = (qubits,)
        self.gates.append(Gate(name, tuple(qubits), param))
        return self

    def extend(self, gates):
        self.gates.extend(gates)
        return self

    def copy(self, gates=None):
        return Circuit(self.n, list(self.gates if gates is None else gates),
                       self.global_phase, None if self.output_map is None else list(self.output_map))

    def __len__(self):
        return len(self.gates)

    def __iter__(self):
        return iter(self.gates)

    def inverse(self):
        from .gates import inverse
        return Circuit(self.n, [inverse(g) for g in reversed(self.gates)], -self.global_phase)

    # -- analysis -----------------------------------------------------------
    def unitary(self) -> np.ndarray:
        n, N = self.n, 2 ** self.n
        t = np.eye(N, dtype=complex).reshape((2,) * n + (N,))
        for g in self.gates:
            t = apply_gate(t, g, n)
        u = t.reshape(N, N) * np.exp(1j * self.global_phase)
        if self.output_map is not None:  # relabel qubits so output bit a sits on qubit a
            perm = np.zeros(N, dtype=int)
            for x in range(N):
                y = 0
                for a, q in enumerate(self.output_map):
                    y |= ((x >> q) & 1) << a
                perm[x] = y
            u2 = np.zeros_like(u)
            u2[perm, :] = u
            u = u2
        return u

    def counts(self) -> Counter:
        return Counter(g.name for g in self.gates)

    def two_qubit_count(self) -> int:
        return sum(1 for g in self.gates if len(g.qubits) >= 2)

    def depth(self) -> int:
        level = [0] * self.n
        for g in self.gates:
            d = max(level[q] for q in g.qubits) + 1
            for q in g.qubits:
                level[q] = d
        return max(level, default=0)

    def summary(self) -> str:
        c = ", ".join(f"{k}:{v}" for k, v in sorted(self.counts().items()))
        return f"{len(self)} gates (2q+: {self.two_qubit_count()}), depth {self.depth()} [{c}]"

    # -- output -------------------------------------------------------------
    def to_qasm(self) -> str:
        lines = ["OPENQASM 3.0;", 'include "stdgates.inc";', f"qubit[{self.n}] q;"]
        if abs(self.global_phase) > 1e-12:
            lines.append(f"gphase({self.global_phase!r});")
        for g in self.gates:
            qs = ", ".join(f"q[{q}]" for q in g.qubits)
            if g.name == "MCP":
                lines.append(f"ctrl({len(g.qubits) - 1}) @ p({g.param!r}) {qs};")
            elif g.name == "CCZ":
                lines.append(f"ctrl(2) @ z {qs};")
            else:
                p = f"({g.param!r})" if g.param is not None else ""
                lines.append(f"{QASM_NAMES[g.name]}{p} {qs};")
        if self.output_map is not None:
            lines.append("// output bit a is on qubit: " + str(self.output_map))
        return "\n".join(lines) + "\n"

    def draw(self) -> str:
        """ASCII diagram, qubit 0 on top."""
        n = self.n
        rows = [[f"q{q}: "] for q in range(n)]
        w0 = max(len(r[0]) for r in rows)
        for r in rows:
            r[0] = r[0].rjust(w0)
        for g in self.gates:
            cells = {}
            ph = phase_of(g)
            if g.name == "CX":
                cells = {g.qubits[0]: "●", g.qubits[1]: "X"}
            elif g.name == "CCX":
                cells = {g.qubits[0]: "●", g.qubits[1]: "●", g.qubits[2]: "X"}
            elif g.name == "SWAP":
                cells = {g.qubits[0]: "x", g.qubits[1]: "x"}
            elif ph is not None and len(g.qubits) > 1:
                qs, th = ph
                cells = {q: "●" for q in qs}
                cells[qs[-1]] = "Z" if g.name in ("CZ", "CCZ") else f"P({fmt_angle(th)})"
            else:
                label = g.name if g.param is None else f"{g.name}({fmt_angle(g.param)})"
                cells = {g.qubits[0]: label}
            width = max(len(c) for c in cells.values()) + 2
            lo, hi = min(g.qubits), max(g.qubits)
            for q in range(n):
                if q in cells:
                    s = cells[q].center(width, "─")
                elif lo < q < hi:
                    s = "┼".center(width, "─")
                else:
                    s = "─" * width
                rows[q].append(s)
        return "\n".join("".join(r) for r in rows)


def phase_key(u: np.ndarray, decimals: int = 7) -> bytes:
    """Hashable fingerprint of a unitary up to global phase."""
    flat = u.reshape(-1)
    idx = int(np.argmax(np.abs(flat) > 1e-6))
    v = u * (abs(flat[idx]) / flat[idx])
    v = np.round(v, decimals) + (0.0 + 0.0j)
    return v.tobytes()


def distance(u: np.ndarray, v: np.ndarray) -> float:
    """Frobenius distance between u and v after optimal global phase alignment."""
    ip = np.vdot(v, u)
    ph = ip / abs(ip) if abs(ip) > 1e-12 else 1.0
    return float(np.linalg.norm(u - ph * v))


def fidelity(u: np.ndarray, v: np.ndarray) -> float:
    """|Tr(U^dagger V)| / N."""
    return float(abs(np.vdot(u, v)) / u.shape[0])


def equivalent(u, v, tol=1e-8) -> bool:
    return distance(u, v) < tol * math.sqrt(u.size)
