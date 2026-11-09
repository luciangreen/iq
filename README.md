# iq

**iq** works out quantum circuits for algorithms, starting from first principles. You give it the *definition* of a transform, not a circuit. iq expands the transform's phase into a Boolean polynomial, cancels every term whose coefficient is a whole number of turns, and reads the circuit off the terms that are left. The result is then simplified by cancelling, merging and commuting gates, optionally lowered to CNOT + phase, and checked against the original definition.

iq is written in SWI-Prolog and uses exact rational arithmetic (`rdiv`), so every cancellation is exact, not a floating-point approximation.

```
swipl prolog/cli.pl qft 4
swipl prolog/cli.pl toffoli --basis cx --qasm
swipl prolog/cli.pl qft 8 --approx 0.05 --quiet
swipl prolog/cli.pl targets
swipl -g run_tests -t halt tests/test_iq.pl
```

Interactive:

```prolog
?- use_module(prolog/iq).
?- iq(qft, 4).
?- iq(toffoli, _, [basis(cx)]).
?- derive(qft, 5, [], derivation(Name, Method, Raw, Circ, Trace, Error)).
```

## The idea: cancelling parts

The QFT is defined by `U|k⟩ = N^-½ Σ_j e^{2πi jk/N} |j⟩`. Measure the phase in turns and write both integers in binary:

```
φ(j,k) = jk/N = Σ_a Σ_b 2^(a+b-n) · j_a·k_b
```

That gives n² terms. Because e^{2πi·m} = 1 for every integer m, each term with a+b ≥ n has an integer coefficient, so it **cancels**. For n = 4 iq reports:

```
φ has 16 terms; reduce mod 1 (e^{2πi·m} = 1 for integer m):
  cancelled 6 integer terms: j1·k3:1, j2·k2:1, j2·k3:2, j3·k1:1, j3·k2:2, j3·k3:4
  surviving 10 terms: φ ≡ 1/16·j0·k0 + 1/8·j0·k1 + 1/4·j0·k2 + 1/2·j0·k3 + …
```

The surviving terms are then read as gates:

| surviving term | meaning | gate |
|---|---|---|
| `½ · j_a·k_b` | (-1)^{j_a k_b}, the Hadamard kernel | `H` on the qubit holding `k_b`, which then holds `j_a` |
| `c · j_a·k_b` (other c) | phase that couples an output bit to an input bit | `CP(2πc)`, applied *after* H on j_a's qubit and *before* H on k_b's qubit |
| `c · k_a·k_b·…` / `c · j_a·j_b·…` | phase on inputs only / outputs only | (multi-)controlled phase before / after the Hadamards |
| leftover pairing a → π(a) | output bits land permuted (bit reversal for the QFT) | SWAP network, or `--no-swaps` to relabel |

The "after / before" rules become ordering constraints, and a topological sort of those constraints gives the gate order. For the QFT this reproduces the textbook circuit: n H gates, n(n−1)/2 controlled phases and ⌊n/2⌋ SWAPs. Nothing about that circuit was put in by hand. The symbolic route never builds a matrix, so `qft 20` works. Verification is only automatic up to n = 6.

## Beyond the QFT: recognising structure in any matrix

For other targets, iq builds the matrix straight from its definition (in `targets.pl`) and tries these steps in order:

1. **Diagonal in some Hadamard frame** `H_S·U·H_S`. A Möbius transform gives the exact multilinear phase polynomial, and every term becomes a phase gate. A *fixed-polarity* search (conjugating by X) picks the form with the fewest terms. This finds CNOT = H·CZ·H, Toffoli = H·CCZ·H, and Grover diffusion = H X·MCZ·X H.
2. **Permutation** (classical reversible). Shortest-circuit search over {X, CX}, then {X, CX, CCX}. This finds SWAP = 3 CX and Fredkin = CX·CCX·CX.
3. **Flat amplitudes** `|U[j,k]| = 1/√N` in some frame. Möbius transform over the 2n bits j, k, then the Hadamard-kernel reading above. This covers H⊗n, the QFT, inverse QFT and similar transforms.
4. **Otherwise**, breadth-first search over Clifford+T+CX, deduplicated up to global phase (finds iSWAP).

## Lowering to CNOT + phase (`--basis cx`)

Each block of consecutive diagonal gates is merged into a single phase polynomial and rewritten in the parity basis:

```
x·y   = ½ (x + y − x⊕y)
x·y·z = ¼ (x + y + z − x⊕y − x⊕z − y⊕z + x⊕y⊕z)
```

The result is reduced mod 1 again, which cancels more terms, and built as a CNOT ladder per parity. The parities are put in Gray-code order so that neighbouring ladders cancel in the simplifier. From the definition of Toffoli, this gives the standard **6-CNOT, 7-T** circuit.

## Simplifier

The simplifier repeats these rewrites until nothing changes:

* cancel `H·H`, `X·X`, `CX·CX`, `SWAP·SWAP`
* drop phases equal to 0 mod 1
* merge phases on the same qubits
* commute gates past each other (disjoint qubits, diagonal gates, CX rules)
* `H·X·H → Z`, `H·Z·H → X`

`--approx ε` also drops rotations smaller than ε radians (the approximate QFT) and reports the fidelity.

## Layout

```
prolog/
  poly.pl      Boolean phase polynomials, mod-1 cancellation, parity basis
  gates.pl     gate terms: h/1 x/1 p(Qubits,Turns) cx/2 ccx/3 swap/2
  synth.pl     reading circuits off reduced polynomials (Hadamard kernels, ordering, SWAPs)
  simplify.pl  peephole rewriting, approximation, CNOT+phase lowering
  derive.pl    recognition pipeline, Möbius transform, symbolic QFT, verification
  search.pl    breadth-first search up to global phase
  sim.pl       state-vector simulation, distance, fidelity
  targets.pl   target unitaries from their definitions
  output.pl    summary, ASCII drawing, OpenQASM 3
  iq.pl        library entry point;  cli.pl  command line
tests/test_iq.pl   plunit suite
reference/iq/      Python implementation of the same algorithm (numpy), used as an oracle
```

Phase gates are stored as `p(Qubits, Turns)` with exact rational turns, so Z, S, T, CZ and CCZ are all `p/2` with turns 1/2, 1/4, 1/8, …. They get their usual names only when printed.

## Status

The Python reference in `reference/` has been run. Every target above checks out exactly against its definition: Toffoli lowers to 6 CX + 7 T, and the approximate QFT on 6 qubits keeps fidelity 0.9895 with 3 rotations removed. The Prolog version implements the same algorithm and comes with a `plunit` suite. It was written without a Prolog system available, so run the tests first:

```
swipl -g run_tests -t halt tests/test_iq.pl
```

Python reference: `cd reference && python3 -c "from iq.derive import derive; d = derive('qft', 4); print(d.report())"`
