"""Peephole simplification: cancel, merge, commute, fuse; plus approximation
and lowering of diagonal blocks to CNOT + phase (parity networks)."""
from __future__ import annotations

import math
from itertools import product
from fractions import Fraction

import numpy as np

from .circuit import Circuit, phase_key
from .gates import (Gate, DAGGER, SELF_INVERSE, TAU, gate_matrix, is_diagonal,
                    norm_angle, phase_gate, phase_of, fmt_angle)
from .phasepoly import Poly, to_parity_basis, mod1, is_zero, nice, mono_str

PRETTY = {1: {math.pi: "Z", math.pi / 2: "S", -math.pi / 2: "Sdg",
              math.pi / 4: "T", -math.pi / 4: "Tdg"},
          2: {math.pi: "CZ"}, 3: {math.pi: "CCZ"}}


def to_phase_form(g: Gate) -> Gate:
    ph = phase_of(g)
    return phase_gate(*ph) if ph is not None else g


def pretty(g: Gate) -> Gate:
    ph = phase_of(g)
    if ph is None:
        return g
    qs, th = ph
    th = norm_angle(th)
    for ang, name in PRETTY.get(len(qs), {}).items():
        if abs(th - ang) < 1e-10:
            return Gate(name, qs)
    return phase_gate(qs, th)


def _zero_angle(th):
    return abs(norm_angle(th)) < 1e-10


def commutes(g: Gate, h: Gate) -> bool:
    if not set(g.qubits) & set(h.qubits):
        return True
    if is_diagonal(g) and is_diagonal(h):
        return True
    for a, b in ((g, h), (h, g)):
        if a.name == "CX" and b.name == "CX":
            return a.qubits[0] != b.qubits[1] and b.qubits[0] != a.qubits[1]
        if a.name == "CX" and is_diagonal(b):
            return a.qubits[1] not in b.qubits
        if a.name == "CX" and b.name == "X" and len(b.qubits) == 1:
            return b.qubits[0] == a.qubits[1]
    return False


def _inverse_pair(g: Gate, h: Gate) -> bool:
    if g.qubits != h.qubits:
        return False
    if g.name in SELF_INVERSE and g.name == h.name:
        return True
    if DAGGER.get(g.name) == h.name:
        return True
    if g.name == h.name == "RZ":
        return _zero_angle(g.param + h.param)
    return False


def cancel_and_merge(gates: list, log: list) -> bool:
    changed = False
    i = 0
    while i < len(gates):
        g = gates[i]
        ph = phase_of(g)
        if (ph is not None and _zero_angle(ph[1])) or g.name == "I":
            log.append(f"drop identity {g}")
            del gates[i]
            changed = True
            continue
        j, restart = i + 1, False
        while j < len(gates):
            h = gates[j]
            if _inverse_pair(g, h):
                log.append(f"cancel {g} · {h}")
                del gates[j]
                del gates[i]
                changed = restart = True
                break
            ph2 = phase_of(h)
            if ph is not None and ph2 is not None and ph[0] == ph2[0]:
                th = ph[1] + ph2[1]
                log.append(f"merge {g} · {h} -> phase {fmt_angle(norm_angle(th))}")
                del gates[j]
                changed = True
                if _zero_angle(th):
                    del gates[i]
                    restart = True
                    break
                g = gates[i] = phase_gate(ph[0], th)
                ph = phase_of(g)
                continue
            if g.name == "RZ" and h.name == "RZ" and g.qubits == h.qubits:
                log.append(f"merge {g} · {h}")
                del gates[j]
                g = gates[i] = Gate("RZ", g.qubits, g.param + h.param)
                changed = True
                continue
            if commutes(g, h):
                j += 1
                continue
            break
        if restart:
            i = max(i - 1, 0)
        else:
            i += 1
    return changed


_TABLE = None


def _single_qubit_table(max_len=3):
    """Shortest word over a Clifford+T alphabet for each 2x2 unitary (up to phase)."""
    global _TABLE
    if _TABLE is None:
        alphabet = ["H", "X", "Y", "Z", "S", "Sdg", "T", "Tdg"]
        table = {phase_key(np.eye(2)): ()}
        for L in range(1, max_len + 1):
            for word in product(alphabet, repeat=L):
                m = np.eye(2)
                for w in word:
                    m = gate_matrix(Gate(w, (0,))) @ m
                table.setdefault(phase_key(m), word)
        _TABLE = table
    return _TABLE


def fuse_single_qubit(gates: list, log: list) -> bool:
    table = _single_qubit_table()
    runs, replace, removed = {}, {}, set()

    def flush(q):
        idxs = runs.pop(q, [])
        if len(idxs) < 2:
            return
        m = np.eye(2)
        for i in idxs:
            m = gate_matrix(Gate(gates[i].name, (0,), gates[i].param)) @ m
        word = table.get(phase_key(m))
        if word is not None and len(word) < len(idxs):
            old = " ".join(gates[i].name for i in idxs)
            log.append(f"fuse q{q}: {old} -> {' '.join(word) or 'I'}")
            replace[idxs[0]] = [Gate(w, (q,)) for w in word]
            removed.update(idxs)

    for i, g in enumerate(gates):
        if len(g.qubits) == 1:
            runs.setdefault(g.qubits[0], []).append(i)
        else:
            for q in g.qubits:
                flush(q)
    for q in list(runs):
        flush(q)
    if not removed:
        return False
    new = []
    for i, g in enumerate(gates):
        if i in replace:
            new.extend(replace[i])
        elif i not in removed:
            new.append(g)
    gates[:] = new
    return True


def simplify(circ: Circuit, log: list | None = None) -> Circuit:
    log = [] if log is None else log
    gates = [to_phase_form(g) for g in circ.gates]
    while True:
        changed = cancel_and_merge(gates, log)
        gates = [pretty(g) for g in gates]
        changed |= fuse_single_qubit(gates, log)
        gates = [to_phase_form(g) for g in gates]
        if not changed:
            break
    return circ.copy([pretty(g) for g in gates])


def approximate(circ: Circuit, eps: float, log: list | None = None) -> Circuit:
    """Drop phase rotations with |theta| < eps (e.g. the approximate QFT)."""
    log = [] if log is None else log
    keep = []
    for g in circ.gates:
        ph = phase_of(g)
        if ph is not None and abs(norm_angle(ph[1])) < eps:
            log.append(f"drop small rotation {g}")
        else:
            keep.append(g)
    return circ.copy(keep)


def _gray_rank(mask: int) -> int:
    b = 0
    while mask:
        b ^= mask
        mask >>= 1
    return b


def parity_network(qubit_phases: Poly, log: list) -> list:
    """Synthesize a diagonal phase polynomial over vars 'x<q>' as CNOT + P."""
    parities = to_parity_basis(qubit_phases)
    kept, ncancel = {}, 0
    for t, c in parities.items():
        if not t:
            continue
        r = mod1(nice(c))
        if is_zero(r):
            ncancel += 1
        else:
            kept[t] = r
    log.append(f"parity basis: {len(parities)} parities, {ncancel} cancel mod 1, {len(kept)} remain: "
               + ", ".join(f"{c}·[{'⊕'.join(sorted(t, key=lambda v: int(v[1:])))}]"
                           for t, c in sorted(kept.items(), key=lambda kv: (len(kv[0]), sorted(kv[0])))))
    order = []
    for t, c in kept.items():
        qs = sorted(int(v[1:]) for v in t)
        tgt = qs[-1]
        mask = sum(1 << q for q in qs[:-1])
        order.append((tgt, _gray_rank(mask), qs, c))
    order.sort()
    gates = []
    for tgt, _, qs, c in order:
        ladder = [Gate("CX", (q, tgt)) for q in qs[:-1]]
        gates += ladder + [phase_gate((tgt,), TAU * float(c))] + ladder[::-1]
    return gates


def lower_to_cx(circ: Circuit, log: list | None = None) -> Circuit:
    """Replace each maximal block of consecutive diagonal phase gates with a
    parity (CNOT + P) network derived from the block's combined phase polynomial."""
    log = [] if log is None else log
    out, block = [], []

    def flush():
        if not block:
            return
        if all(len(g.qubits) == 1 for g in block):
            out.extend(block)
        else:
            poly = Poly()
            for g in block:
                qs, th = phase_of(g)
                mono = Poly.const(nice(th / TAU))
                for q in qs:
                    mono = mono * Poly.var(f"x{q}")
                poly = poly + mono
            log.append("diagonal block " + " ".join(str(g) for g in block) + f"  =>  φ = {poly}")
            out.extend(parity_network(poly, log))
        block.clear()

    for g in circ.gates:
        if phase_of(g) is not None:
            block.append(g)
        else:
            flush()
            out.append(g)
    flush()
    return circ.copy(out)
