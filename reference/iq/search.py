"""Bidirectional (meet-in-the-middle) breadth-first search for a gate sequence
equal to a target unitary up to global phase.  Used when no algebraic
structure is recognised, and to confirm minimal forms of small circuits."""
from __future__ import annotations

from itertools import combinations, permutations

import numpy as np

from .circuit import Circuit, apply_gate, phase_key
from .gates import ARITY, SYMMETRIC, Gate, inverse, parse_gate_spec

DEFAULT_LIBRARY = ["H", "S", "Sdg", "T", "Tdg", "X", "CX"]


def instantiate(library, n):
    gates = []
    for spec in library:
        name, param = parse_gate_spec(spec) if isinstance(spec, str) else spec
        k = ARITY.get(name, 2 if name == "CP" else 1)
        if k > n:
            continue
        tuples = combinations(range(n), k) if name in SYMMETRIC else permutations(range(n), k)
        for qs in tuples:
            gates.append(Gate(name, qs, param))
    return gates


def _apply(u, g, n):
    N = 2 ** n
    return apply_gate(u.reshape((2,) * n + (N,)), g, n).reshape(N, N)


def _bfs(start, gates, n, depth, max_nodes):
    """layers: key -> (unitary, path)."""
    seen = {phase_key(start): (start, [])}
    frontier = [(start, [])]
    for _ in range(depth):
        nxt = []
        for u, path in frontier:
            for g in gates:
                if path and inverse(g) == path[-1]:
                    continue
                v = _apply(u, g, n)
                key = phase_key(v)
                if key not in seen:
                    seen[key] = (v, path + [g])
                    nxt.append((v, path + [g]))
                    if len(seen) >= max_nodes:
                        return seen
        frontier = nxt
    return seen


def search(target: np.ndarray, n: int, library=None, max_depth=8, max_nodes=300_000):
    """Return the shortest Circuit (<= max_depth gates) implementing target, or None."""
    gates = instantiate(library or DEFAULT_LIBRARY, n)
    N = 2 ** n
    tkey = phase_key(target)
    # iterative deepening of the two halves keeps the answer minimal
    for d in range(0, max_depth + 1):
        fd, bd = (d + 1) // 2, d // 2
        fwd = _bfs(np.eye(N, dtype=complex), gates, n, fd, max_nodes)
        if tkey in fwd and len(fwd[tkey][1]) <= d:
            return Circuit(n, fwd[tkey][1])
        # backward: B = g_m^dag ... g_1^dag T ; if B == F then T = g_1 .. g_m F
        inv = [inverse(g) for g in gates]
        bwd = _bfs(target.astype(complex), inv, n, bd, max_nodes)
        best = None
        for key, (_, bpath) in bwd.items():
            if key in fwd:
                fpath = fwd[key][1]
                full = fpath + [inverse(g) for g in reversed(bpath)]
                if best is None or len(full) < len(best):
                    best = full
        if best is not None and len(best) <= d:
            return Circuit(n, best)
    return None
