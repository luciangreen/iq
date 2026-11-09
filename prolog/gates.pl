:- encoding(utf8).
:- module(gates, [gate_qubits/2, mkphase/3, is_phase/1, gate_label/2, gate_string/2,
                  angle_string/2, inverse_gate/2, qasm_line/2]).
/** <module> Gate terms

Conventions: basis index k = sum_b k_b * 2^b, qubit b holds bit k_b.

  h(Q)            Hadamard
  x(Q)            Pauli X
  p(Qs, T)        phase e^{2 pi i T} when ALL qubits in Qs are 1 (Qs sorted);
                  T in turns (exact rational).  P, CP, MCP and the named
                  Z = p([q],1/2), S = 1/4, T = 1/8, CZ = p([a,b],1/2),
                  CCZ = p([a,b,c],1/2) are all p/2 internally.
  cx(C, T)        CNOT
  ccx(A, B, T)    Toffoli, A < B
  swap(A, B)      SWAP, A < B
*/
:- use_module(library(lists)).
:- use_module(poly).

gate_qubits(h(Q), [Q]).
gate_qubits(x(Q), [Q]).
gate_qubits(p(Qs, _), Qs).
gate_qubits(cx(C, T), [C, T]).
gate_qubits(ccx(A, B, T), [A, B, T]).
gate_qubits(swap(A, B), [A, B]).

is_phase(p(_, _)).

mkphase(Qs0, T0, p(Qs, T)) :-
    msort(Qs0, Qs),
    snap(T0, T1),
    mod1(T1, T).

inverse_gate(p(Qs, T), p(Qs, T1)) :- !, T1 is -T.
inverse_gate(G, G).

%!  gate_label(+Gate, -Label)
%   Display name, recognising the named phase gates.
gate_label(p([_], T), z)   :- T =:= 1 rdiv 2, !.
gate_label(p([_], T), s)   :- T =:= 1 rdiv 4, !.
gate_label(p([_], T), sdg) :- T =:= -1 rdiv 4, !.
gate_label(p([_], T), t)   :- T =:= 1 rdiv 8, !.
gate_label(p([_], T), tdg) :- T =:= -1 rdiv 8, !.
gate_label(p([_, _], T), cz) :- T =:= 1 rdiv 2, !.
gate_label(p([_, _, _], T), ccz) :- T =:= 1 rdiv 2, !.
gate_label(p([_], _), p) :- !.
gate_label(p([_, _], _), cp) :- !.
gate_label(p(_, _), mcp) :- !.
gate_label(G, L) :- functor(G, L, _).

%!  angle_string(+Turns, -S)   e.g. 1/8 turn -> "pi/4"
angle_string(T, S) :-
    (   rational(T)
    ->  F is 2*T,
        rational(F, N, D),
        (   N =:= 0 -> S = "0"
        ;   (   N =:= 1 -> H = "pi" ; N =:= -1 -> H = "-pi" ; format(string(H), "~dpi", [N]) ),
            (   D =:= 1 -> S = H ; format(string(S), "~w/~d", [H, D]) )
        )
    ;   R is 2*pi*T,
        format(string(S), "~4f", [R])
    ).

gate_string(G, S) :-
    gate_label(G, L),
    gate_qubits(G, Qs),
    findall(QS, ( member(Q, Qs), format(string(QS), " q~d", [Q]) ), QSs),
    atomic_list_concat(QSs, QA),
    (   G = p(_, T), memberchk(L, [p, cp, mcp])
    ->  angle_string(T, AS),
        format(string(S), "~w(~w)~w", [L, AS, QA])
    ;   format(string(S), "~w~w", [L, QA])
    ).

%!  qasm_line(+Gate, -Line)   OpenQASM 3
qasm_line(G, Line) :-
    gate_label(G, L),
    gate_qubits(G, Qs),
    findall(QS, ( member(Q, Qs), format(string(QS), "q[~d]", [Q]) ), QSs),
    atomic_list_concat(QSs, ', ', QA),
    (   L == ccz
    ->  format(string(Line), "ctrl(2) @ z ~w;", [QA])
    ;   G = p(_, T), memberchk(L, [p, cp, mcp])
    ->  R is 2*pi*T,
        (   L == mcp
        ->  length(Qs, K), K1 is K - 1,
            format(string(Line), "ctrl(~d) @ p(~15g) ~w;", [K1, R, QA])
        ;   format(string(Line), "~w(~15g) ~w;", [L, R, QA])
        )
    ;   format(string(Line), "~w ~w;", [L, QA])
    ).
