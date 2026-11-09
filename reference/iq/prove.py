"""Exact equivalence proofs by sum-over-paths (path sums).

A circuit on n qubits acts as
    U|x> = 2^{-h/2} sum_{y in {0,1}^h} e^{2 pi i phi(x,y)} |f(x,y)>
with phi a Boolean phase polynomial (turns, mod 1) and f a vector of GF(2)
polynomials.  To prove C == U we build the path sum of C, append the path
sum of U^dagger (written from U's *definition*), and reduce with two sound
rewrite rules until the identity remains:

  [Elim] a path variable y appearing nowhere: sum_y 1 = 2       (h -= 2)
  [HH]   y not in outputs, phi = phi' + 1/2 * y * (v xor Q), v another path
         variable not in Q:  sum_y (-1)^{y (v xor Q)} = 2 [v = Q]
         => drop y, substitute v := Q everywhere              (h -= 2)

Identity = no path variables, outputs f_q = x_q, phase constant (global).
"""
from __future__ import annotations

from fractions import Fraction

from .phasepoly import Poly, mod1, is_zero, mono_str, var_key

HALF = Fraction(1, 2)

# ---- GF(2) polynomials: frozenset of monomials (frozenset of var names) ----
ONE = frozenset([frozenset()])
ZERO = frozenset()


def gvar(v):
    return frozenset([frozenset([v])])


def gxor(a, b):
    return a ^ b


def gand(a, b):
    out = set()
    for m1 in a:
        for m2 in b:
            out ^= {m1 | m2}
    return frozenset(out)


def gsubst(p, v, q):
    out = set()
    for m in p:
        if v in m:
            for t in gand(frozenset([m - {v}]), q):
                out ^= {t}
        else:
            out ^= {m}
    return frozenset(out)


def gstr(p):
    if not p:
        return "0"
    return " ⊕ ".join(mono_str(m) for m in sorted(p, key=lambda m: (len(m), sorted(map(var_key, m)))))


def lift(p) -> Poly:
    """Arithmetic value of a GF(2) polynomial: a xor b = a + b - 2ab."""
    acc = Poly()
    for m in p:
        mp = Poly({m: Fraction(1)})
        acc = acc + mp - acc * mp * 2
    return acc


def psubst(phi: Poly, v, q) -> Poly:
    lq = lift(q)
    out = Poly()
    for m, c in phi.terms.items():
        if v in m:
            out = out + Poly({m - {v}: c}) * lq
        else:
            out = out + Poly({m: c})
    return out


class PathSum:
    def __init__(self, n):
        self.n = n
        self.out = [gvar(f"x{q}") for q in range(n)]
        self.phase = Poly()
        self.h = 0
        self.paths = []
        self._next = 0

    def fresh(self):
        v = f"y{self._next}"
        self._next += 1
        self.paths.append(v)
        self.h += 1
        return v

    def add_phase(self, p: Poly):
        self.phase, _, _ = (self.phase + p).reduce_mod1()

    # gates ---------------------------------------------------------------
    def apply(self, g):
        n, q = g.name, g.qubits
        from .gates import phase_of
        ph = phase_of(g)
        if n == "H":
            y = self.fresh()
            self.add_phase(lift(self.out[q[0]]) * Poly.var(y) * HALF)
            self.out[q[0]] = gvar(y)
        elif n == "X":
            self.out[q[0]] = gxor(self.out[q[0]], ONE)
        elif n == "CX":
            self.out[q[1]] = gxor(self.out[q[1]], self.out[q[0]])
        elif n == "CCX":
            self.out[q[2]] = gxor(self.out[q[2]], gand(self.out[q[0]], self.out[q[1]]))
        elif n == "SWAP":
            a, b = q
            self.out[a], self.out[b] = self.out[b], self.out[a]
        elif ph is not None:
            qs, th = ph
            t = Fraction(th / (2 * 3.141592653589793)).limit_denominator(1 << 24)
            p = Poly.const(t)
            for qq in qs:
                p = p * lift(self.out[qq])
            self.add_phase(p)
        else:
            raise ValueError(f"no path-sum rule for {g}")

    def state(self):
        outs = ", ".join(f"q{q}:{gstr(o)}" for q, o in enumerate(self.out))
        return f"2^(-{self.h}/2) Σ_{{{','.join(self.paths)}}} e^{{2πi[{self.phase}]}} |{outs}⟩"

    # reduction -------------------------------------------------------------
    def in_outputs(self, v):
        return any(v in m for o in self.out for m in o)

    def reduce(self, log):
        changed = True
        while changed:
            changed = False
            for y in list(self.paths):
                if self.in_outputs(y):
                    continue
                ty = {m: c for m, c in self.phase.terms.items() if y in m}
                if not ty:
                    self.paths.remove(y)
                    self.h -= 2
                    log.append(f"[Elim] {y} occurs nowhere: Σ_{y} 1 = 2")
                    changed = True
                    break
                if not all(is_zero(mod1(c - HALF)) for c in ty.values()):
                    continue
                r = frozenset()
                for m in ty:
                    r ^= {m - {y}}
                if ONE <= r and len(r) == 1:
                    raise ValueError("path sum vanishes (not unitary)")
                pick = None
                for v in self.paths:
                    if v != y and frozenset([v]) in r and not any(v in m for m in r if m != frozenset([v])):
                        pick = v
                        break
                if pick is None:
                    continue
                q = r ^ {frozenset([pick])}
                for m in ty:
                    del self.phase.terms[m]
                self.phase = psubst(self.phase, pick, q)
                self.phase, _, _ = self.phase.reduce_mod1()
                self.out = [gsubst(o, pick, q) for o in self.out]
                self.paths.remove(y)
                self.paths.remove(pick)
                self.h -= 2
                log.append(f"[HH] Σ_{y} (-1)^{{{y}·({gstr(r)})}} forces {pick} := {gstr(q) if q else '0'}")
                changed = True
                break

    def is_identity(self):
        const = [m for m in self.phase.terms if m]
        return (not self.paths and self.h == 0 and not const and
                all(o == gvar(f"x{q}") for q, o in enumerate(self.out)))


# ---- inverse of each target, written from its definition -----------------
def _register(ps, outs):
    """sum_b 2^b * lift(out_b)."""
    acc = Poly()
    for b, o in enumerate(outs):
        acc = acc + lift(o) * (2 ** b)
    return acc


def apply_definition_inverse(ps: PathSum, name: str):
    n, N = ps.n, 2 ** ps.n
    o = ps.out
    if name in ("qft", "iqft", "hadamard"):
        K = _register(ps, o) if name != "hadamard" else None
        zs = [ps.fresh() for _ in range(n)]
        if name == "hadamard":
            phi = Poly()
            for a in range(n):
                phi = phi + Poly.var(zs[a]) * lift(o[a]) * HALF
            ps.add_phase(phi * -1)
            desc = "H^n†|k⟩ = 2^(-n/2) Σ_z (-1)^{z·k} |z⟩"
        else:
            J = Poly()
            for a in range(n):
                J = J + Poly.var(zs[a]) * (2 ** a)
            sgn = Fraction(-1 if name == "qft" else 1, N)
            ps.add_phase(J * K * sgn)
            desc = f"{name}†|k⟩ = 2^(-n/2) Σ_z e^{{{'-' if name == 'qft' else '+'}2πi z·k/{N}}} |z⟩"
        ps.out = [gvar(z) for z in zs]
        return desc
    if name in ("cz", "ccz", "mcz"):
        p = Poly.const(-HALF)
        for x in o:
            p = p * lift(x)
        ps.add_phase(p)
        return "phase −½·Π k_b"
    if name == "chirp":
        K = _register(ps, o)
        ps.add_phase(K * K * Fraction(-1, 2 * N))
        return f"phase −k²/{2 * N}"
    if name == "cx":
        ps.out[1] = gxor(o[1], o[0])
        return "|a,b⟩ → |a, b⊕a⟩"
    if name == "toffoli":
        ps.out[2] = gxor(o[2], gand(o[0], o[1]))
        return "|a,b,c⟩ → |a, b, c⊕ab⟩"
    if name == "swap":
        ps.out[0], ps.out[1] = o[1], o[0]
        return "|a,b⟩ → |b,a⟩"
    if name == "fredkin":
        d = gand(o[0], gxor(o[1], o[2]))
        ps.out[1], ps.out[2] = gxor(o[1], d), gxor(o[2], d)
        return "|c,a,b⟩ → |c, a⊕c(a⊕b), b⊕c(a⊕b)⟩"
    if name == "iswap":
        ps.add_phase(lift(gxor(o[0], o[1])) * Fraction(-1, 4))
        ps.out[0], ps.out[1] = ps.out[1], ps.out[0]
        return "|a,b⟩ → i^{-(a⊕b)} |b,a⟩"
    if name == "diffusion":
        from .gates import Gate
        for q in range(n):
            ps.apply(Gate("H", (q,)))
        p = Poly.const(1)
        for x in ps.out:
            p = p * (Poly.const(1) - lift(x))
        ps.add_phase((Poly.const(1) - p) * HALF)
        for q in range(n):
            ps.apply(Gate("H", (q,)))
        return "2|s⟩⟨s|−I = H^n (2|0⟩⟨0|−I) H^n  (self-inverse)"
    raise KeyError(name)


def prove(circuit, name: str):
    """Return (proved: bool, log lines)."""
    log = []
    ps = PathSum(circuit.n)
    for g in circuit.gates:
        ps.apply(g)
    if circuit.output_map is not None:
        ps.out = [ps.out[q] for q in circuit.output_map]
    log.append(f"circuit path sum: {ps.state()}")
    desc = apply_definition_inverse(ps, name)
    log.append(f"compose with U† from the definition: {desc}")
    log.append(f"  {ps.state()}")
    ps.reduce(log)
    ok = ps.is_identity()
    log.append(f"normal form: {ps.state()}")
    gp = ps.phase.terms.get(frozenset(), 0)
    if ok:
        log.append(f"∎ identity: circuit = U exactly (global phase {gp} turn)")
    else:
        log.append("not proved: residual above is what differs (or rules incomplete)")
    return ok, log
