"""First-principles derivation of circuits.

Pipeline for a target U:

1.  Recognise structure.
    * diagonal            U = diag(e^{2πi φ(k)})
    * flat ("Hadamard-type") |U[j,k]| = 1/√N, U[j,k] = e^{2πi φ(j,k)}/√N
    * either of the above after conjugating some qubits by H (e.g. CNOT, Toffoli,
      Grover diffusion)
2.  Expand φ as a multilinear Boolean polynomial (symbolically from a formula,
    or by Möbius transform of the matrix), reduce mod 1 => integer terms cancel.
3.  Read the circuit off the surviving terms:
    * j_a·k_b with coefficient 1/2 is a Hadamard kernel (-1)^{j_a k_b}
    * every other term is a (multi-)controlled phase between the qubits that
      hold its variables at the right moment, which fixes the gate order
    * the pairing j_a <-> k_b leaves a qubit permutation, realised by SWAPs
4.  Simplify (cancellation, merging, commutation, fusion), optionally lower
    to CNOT + phase, optionally approximate, and verify against the definition.
5.  If nothing is recognised: meet-in-the-middle search over a gate library.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction
from itertools import combinations

import numpy as np

from .circuit import Circuit, distance, fidelity
from .gates import TAU, Gate, phase_gate
from .phasepoly import (HALF, Poly, fmt_coef, is_zero, mod1, mono_str, poly_of_diagonal,
                        poly_of_transform, var_key)
from . import simplify as simp
from .search import search, DEFAULT_LIBRARY


class SynthesisError(Exception):
    pass


@dataclass
class Derivation:
    name: str
    circuit: Circuit
    raw: Circuit
    method: str
    trace: list = field(default_factory=list)
    poly: Poly | None = None
    error: float | None = None
    fidelity: float | None = None

    def report(self) -> str:
        return "\n".join(self.trace)


def _split(m):
    J = sorted(int(v[1:]) for v in m if v[0] == "j")
    K = sorted(int(v[1:]) for v in m if v[0] == "k")
    return J, K


# --------------------------------------------------------------------------
# synthesis from a reduced phase polynomial
# --------------------------------------------------------------------------
def synth_diagonal(poly: Poly, n: int, trace: list) -> Circuit:
    c = Circuit(n)
    for m, coef in poly.sorted_terms():
        if not m:
            c.global_phase += TAU * float(coef)
            continue
        qs = [int(v[1:]) for v in m]
        c.gates.append(phase_gate(qs, TAU * float(coef)))
    trace.append(f"  each surviving term c·k_a·k_b… is a phase e^{{2πi c}} on qubits a,b,… "
                 f"=> {len(c)} phase gates")
    return c


def synth_flat(poly: Poly, n: int, trace: list, swaps: bool = True) -> Circuit:
    terms = dict(poly.terms)
    gp = TAU * float(terms.pop(frozenset(), 0))
    # 1. Hadamard kernels: coefficient exactly 1/2 on j_a k_b
    pi = {}
    for m, c in list(terms.items()):
        J, K = _split(m)
        if len(J) == 1 and len(K) == 1 and is_zero(mod1(c - HALF)):
            if J[0] in pi or K[0] in pi.values():
                raise SynthesisError("Hadamard kernels do not form a permutation "
                                     "(needs a CNOT pre-network)")
            pi[J[0]] = K[0]
            del terms[m]
    if len(pi) != n:
        raise SynthesisError(f"only {len(pi)} of {n} output bits have a Hadamard kernel")
    trace.append("  coefficient ½ on j_a·k_b  =>  (-1)^{j_a k_b}: a Hadamard on the qubit "
                 "holding k_b turns it into j_a:")
    trace.append("    " + ", ".join(f"j{a}←H(q{b})" for a, b in sorted(pi.items())))
    # 2. precedence: a term with outputs J and inputs K needs H(pi[a]) done, H(b) not yet
    edges = {q: set() for q in range(n)}
    for m in terms:
        J, K = _split(m)
        for a in J:
            for b in K:
                if pi[a] == b:
                    raise SynthesisError(f"term {mono_str(m)} needs qubit q{b} to hold "
                                         "j and k at once")
                edges[pi[a]].add(b)
    indeg = {q: 0 for q in range(n)}
    for q in edges:
        for b in edges[q]:
            indeg[b] += 1
    order, ready = [], [q for q in range(n) if indeg[q] == 0]
    while ready:
        ready.sort()
        q = ready.pop()  # deterministic: highest available qubit first
        order.append(q)
        for b in edges[q]:
            indeg[b] -= 1
            if indeg[b] == 0:
                ready.append(b)
    if len(order) != n:
        raise SynthesisError("cyclic ordering constraints between controlled phases")
    pos = {q: i for i, q in enumerate(order)}
    trace.append("  a term j_a·k_b must act after H(q{π(a)}) and before H(q_b) "
                 f"=> order {' → '.join(f'H(q{q})' for q in order)}")
    # 3. schedule
    before, after = [], {q: [] for q in range(n)}
    for m, c in terms.items():
        J, K = _split(m)
        if not J:
            before.append((m, c))
        else:
            after[max((pi[a] for a in J), key=pos.get)].append((m, c))
    circ = Circuit(n, global_phase=gp)

    def emit(m, c):
        J, K = _split(m)
        qs = [pi[a] for a in J] + K
        circ.gates.append(phase_gate(qs, TAU * float(c)))

    for m, c in sorted(before, key=lambda t: sorted(map(var_key, t[0]))):
        emit(m, c)
    for q in order:
        circ.add("H", q)
        for m, c in sorted(after[q], key=lambda t: sorted(map(var_key, t[0]))):
            emit(m, c)
    ncp = len(circ) - n
    trace.append(f"  remaining {len(terms)} terms become controlled phases => "
                 f"{n} H + {ncp} phase gates")
    # 4. output permutation
    where = [pi[a] for a in range(n)]
    if all(where[a] == a for a in range(n)):
        trace.append("  output bits already in place: no SWAPs")
    elif swaps:
        nsw = 0
        for a in range(n):
            if where[a] != a:
                b = where.index(a)
                circ.add("SWAP", (a, where[a]))
                where[a], where[b] = a, where[a]
                nsw += 1
        trace.append(f"  output j_a sits on q{{π(a)}} {[pi[a] for a in range(n)]} => {nsw} SWAP(s)")
    else:
        circ.output_map = where
        trace.append(f"  output bits left permuted (no SWAPs): j_a on qubit {where}")
    return circ


# --------------------------------------------------------------------------
# polynomial reduction with a readable account of what cancelled
# --------------------------------------------------------------------------
def reduce_with_trace(poly: Poly, trace: list, label="φ"):
    red, cancelled, kept = poly.reduce_mod1()
    nonconst = [t for t in cancelled if t[0]]
    trace.append(f"  {label} has {len(poly.terms)} terms; reduce mod 1 (e^{{2πi·m}} = 1 for integer m):")
    if nonconst:
        shown = ", ".join(f"{mono_str(m)}:{fmt_coef(c)}" for m, c in
                          sorted(nonconst, key=lambda t: sorted(map(var_key, t[0])))[:24])
        more = "" if len(nonconst) <= 24 else f", … (+{len(nonconst) - 24})"
        trace.append(f"    cancelled {len(nonconst)} integer terms: {shown}{more}")
    else:
        trace.append("    nothing cancels")
    trace.append(f"    surviving {len([k for k in kept if k[0]])} terms: φ ≡ {red}")
    return red


# --------------------------------------------------------------------------
# post-processing shared by every route
# --------------------------------------------------------------------------
def finish(name, raw, method, trace, target=None, poly=None, simplify=True,
           basis="native", approx=None):
    circ = raw
    trace.append(f"Raw circuit: {raw.summary()}")
    if basis == "cx":
        log = []
        circ = simp.lower_to_cx(circ, log)
        trace.append("Lower diagonal blocks to CNOT + phase via the parity basis "
                     "x·y = ½(x + y − x⊕y):")
        trace += ["  " + s for s in log]
    if simplify:
        log = []
        before = len(circ)
        circ = simp.simplify(circ, log)
        trace.append(f"Simplify: {before} → {len(circ)} gates")
        trace += ["  " + s for s in log[:40]]
        if len(log) > 40:
            trace.append(f"  … {len(log) - 40} more rewrites")
    if approx:
        log = []
        circ = simp.approximate(circ, approx, log)
        trace.append(f"Approximate (drop rotations < {approx:g} rad): {len(log)} gates removed")
    d = Derivation(name, circ, raw, method, trace, poly)
    if target is not None:
        u = circ.unitary()
        d.error = distance(u, target)
        d.fidelity = fidelity(target, u)
        ok = "exact" if d.error < 1e-7 else "approximate"
        trace.append(f"Verify against definition: ‖U_circ − e^{{iγ}}U‖ = {d.error:.2e}, "
                     f"fidelity {d.fidelity:.6f} ({ok})")
    trace.append(f"Final circuit: {circ.summary()}")
    return d


# --------------------------------------------------------------------------
# symbolic route: the QFT from its defining formula (scales to large n)
# --------------------------------------------------------------------------
def derive_qft(n: int, inverse=False, shift=0, swaps=True, verify=None, **kw) -> Derivation:
    N = 2 ** n
    sgn = -1 if inverse else 1
    trace = [f"{'Inverse ' if inverse else ''}QFT on {n} qubits from its definition:",
             f"  U|k⟩ = N^-½ Σ_j e^{{2πi·{'-' if inverse else ''}j(k{'+' + str(shift) if shift else ''})/N}} |j⟩,"
             f"  N = {N}",
             "  phase in turns φ(j,k) = " + ("-" if inverse else "") +
             f"j·(k{'+' + str(shift) if shift else ''})/{N};  write j = Σ 2^a j_a, k = Σ 2^b k_b"]
    J, K = Poly.register("j", n), Poly.register("k", n)
    phi = J * (K + shift) * Fraction(sgn, N)
    red = reduce_with_trace(phi, trace)
    trace.append("Read off the circuit:")
    raw = synth_flat(red, n, trace, swaps=swaps)
    from .targets import qft as qft_def
    target = None
    if verify if verify is not None else n <= 9:
        target = qft_def(n, shift)
        target = target.conj().T if inverse else target
        if inverse and shift:
            target = None
    return finish(f"{'i' if inverse else ''}qft{n}", raw, "symbolic phase polynomial",
                  trace, target, red, **kw)


# --------------------------------------------------------------------------
# numeric route: any matrix
# --------------------------------------------------------------------------
def _h_frame(n, S):
    c = Circuit(n)
    for q in S:
        c.add("H", q)
    return c


def _is_diag(u, tol=1e-9):
    return np.allclose(u, np.diag(np.diag(u)), atol=tol)


def _is_flat(u, n, tol=1e-9):
    return np.allclose(np.abs(u), 1 / math.sqrt(2 ** n), atol=tol)


def _is_perm(u):
    a = np.abs(u)
    return bool(np.allclose(a, np.round(a)) and np.allclose(u.imag, 0) and np.allclose(u.real, a))


def _best_polarity(d, n, trace):
    """Fixed-polarity Reed–Muller: choose X-conjugation minimising the term count."""
    best = None
    for r in range(n + 1):
        for S in combinations(range(n), r):
            mask = sum(1 << q for q in S)
            dd = np.array([d[x ^ mask] for x in range(2 ** n)])
            red, _, _ = poly_of_diagonal(dd, n).reduce_mod1()
            nterms = len([m for m in red.terms if m])
            if best is None or nterms < best[0]:
                best = (nterms, S, dd)
    return best


def derive_matrix(u: np.ndarray, name="U", swaps=True, search_depth=8, library=None,
                  **kw) -> Derivation:
    u = np.asarray(u, dtype=complex)
    N = u.shape[0]
    n = int(round(math.log2(N)))
    if 2 ** n != N or u.shape != (N, N):
        raise ValueError("matrix must be 2^n x 2^n")
    if not np.allclose(u.conj().T @ u, np.eye(N), atol=1e-8):
        raise ValueError("matrix is not unitary")
    subsets = [S for r in range(n + 1) for S in combinations(range(n), r)]
    frames = []
    for S in subsets:
        hf = _h_frame(n, S).unitary()
        frames.append((S, hf @ u @ hf))
    # 1. diagonal in some Hadamard frame
    for S, v in frames:
        if _is_diag(v):
            frame = f" in the Hadamard frame on {['q%d' % q for q in S]} (U = H_S·V·H_S)" if S else ""
            trace = [f"{name}: diagonal{frame}: V = diag(e^{{2πi φ(k)}})"]
            d = np.diag(v)
            nterms, X, dd = _best_polarity(d, n, trace)
            if X:
                trace.append(f"  fixed-polarity choice: conjugate by X on {['q%d' % q for q in X]} "
                             f"(literals 1−k_b) => fewest terms ({nterms})")
            poly = poly_of_diagonal(dd, n)
            trace.append(f"  Möbius transform over {n} Boolean variables: φ = {poly}")
            red = reduce_with_trace(poly, trace)
            core = synth_diagonal(red, n, trace)
            pre = [Gate("H", (q,)) for q in S] + [Gate("X", (q,)) for q in X]
            raw = Circuit(n, pre + core.gates + pre[::-1], core.global_phase)
            return finish(name, raw, "diagonal phase polynomial" + (" + H frame" if S else ""),
                          trace, u, red, **kw)
    # 2. classical reversible (permutation) targets: reversible-gate search
    if _is_perm(u):
        for lib in (["X", "CX"], ["X", "CX", "CCX"]):
            found = search(u, n, lib, max_depth=6)
            if found is not None:
                trace = [f"{name}: a permutation of basis states (classical reversible); "
                         f"shortest circuit over {lib}: {len(found)} gates"]
                return finish(name, found, "reversible search", trace, u, None, **kw)
    # 3. Hadamard-type in some frame
    for S, v in frames:
        if _is_flat(v, n):
            frame = f" in the Hadamard frame on {['q%d' % q for q in S]} (U = H_S·V·H_S)" if S else ""
            trace = [f"{name}: all amplitudes 1/√{N}{frame}: V[j,k] = e^{{2πi φ(j,k)}}/√N"]
            poly = poly_of_transform(v, n)
            trace.append(f"  Möbius transform over {2 * n} Boolean variables (j and k bits)")
            red = reduce_with_trace(poly, trace)
            trace.append("Read off the circuit:")
            try:
                core = synth_flat(red, n, trace, swaps=swaps)
            except SynthesisError as e:
                continue
            pre = [Gate("H", (q,)) for q in S]
            raw = Circuit(n, pre + core.gates + pre[::-1], core.global_phase, core.output_map)
            return finish(name, raw, "Hadamard-type phase polynomial" + (" + H frame" if S else ""),
                          trace, u, red, **kw)
    trace = [f"{name}: no diagonal or Hadamard-type structure in any H frame; "
             f"meet-in-the-middle search (≤{search_depth} gates)"]
    found = search(u, n, library or (DEFAULT_LIBRARY + ['CCX']), max_depth=search_depth)
    if found is None:
        raise SynthesisError(f"no circuit of ≤{search_depth} gates found for {name}")
    trace.append(f"  found {len(found)}-gate circuit")
    return finish(name, found, "search", trace, u, None, **kw)


def derive(name: str, n=None, **kw) -> Derivation:
    from .targets import build
    if name in ("qft", "iqft"):
        return derive_qft(n or 3, inverse=(name == "iqft"), **kw)
    u, n = build(name, n)
    return derive_matrix(u, f"{name}{n if name in ('hadamard','chirp','diffusion','mcz') else ''}", **kw)
