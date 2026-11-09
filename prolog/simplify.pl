:- encoding(utf8).
:- module(simplify, [simplify_gates/3, approximate/4, lower_to_cx/3, commutes/2]).
/** <module> Simplification

Rewrites applied to a fixpoint (each strictly shortens the circuit):

  * drop identity        p(Qs, T) with T = 0 mod 1
  * cancel               G ... G^-1 -> nothing   (H H, X X, CX CX, SWAP SWAP)
  * merge                p(Qs,a) ... p(Qs,b) -> p(Qs, a+b)   (phases add)
  * commute past         gates on disjoint qubits, diagonal gates, CX rules
  * conjugate            H X H -> Z,  H Z H -> X

Lowering: every maximal block of diagonal phase gates is rewritten as one
phase polynomial, converted to the parity (XOR) basis with
x*y = (x + y - x xor y)/2, reduced mod 1 (more cancellation) and
synthesised as a CNOT + P parity network, Gray-code ordered so neighbouring
CNOT ladders cancel.
*/
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(pairs)).
:- use_module(poly).
:- use_module(gates).

%!  commutes(+G, +H) is semidet.
commutes(G, H) :-
    gate_qubits(G, A), gate_qubits(H, B),
    \+ ( member(X, A), memberchk(X, B) ), !.
commutes(p(_, _), p(_, _)) :- !.
commutes(cx(C1, T1), cx(C2, T2)) :- !, C1 \== T2, C2 \== T1.
commutes(cx(_, T), p(Qs, _)) :- !, \+ memberchk(T, Qs).
commutes(p(Qs, _), cx(_, T)) :- !, \+ memberchk(T, Qs).
commutes(cx(_, T), x(Q)) :- !, Q == T.
commutes(x(Q), cx(_, T)) :- !, Q == T.

%!  combine(+G, +H, -Result) is semidet.   G...H -> Result (list)
combine(h(Q), h(Q), []).
combine(x(Q), x(Q), []).
combine(cx(C, T), cx(C, T), []).
combine(ccx(A, B, T), ccx(A, B, T), []).
combine(swap(A, B), swap(A, B), []).
combine(p(Qs, T1), p(Qs, T2), R) :-
    T is T1 + T2, mod1(T, T3),
    (   is_zero(T3) -> R = [] ; R = [p(Qs, T3)] ).

partner(G, [H|T], [], H, T) :- combine(G, H, _), !.
partner(G, [H|T], [H|M], P, A) :- commutes(G, H), partner(G, T, M, P, A).

hconj(Q, x(Q), p([Q], Half)) :- Half is 1 rdiv 2.
hconj(Q, p([Q], T), x(Q)) :- T =:= 1 rdiv 2.

rewrite(p(Qs, T), Rest, Rest, Msg) :-
    mod1(T, R), is_zero(R), !,
    gate_string(p(Qs, T), S), format(string(Msg), "drop identity ~w", [S]).
rewrite(G, Rest, Rest2, Msg) :-
    partner(G, Rest, Mid, H, After),
    combine(G, H, R), !,
    append(Mid, R, A1), append(A1, After, Rest2),
    gate_string(G, SG), gate_string(H, SH),
    (   R == []
    ->  format(string(Msg), "cancel ~w · ~w", [SG, SH])
    ;   R = [M], gate_string(M, SM),
        format(string(Msg), "merge ~w · ~w -> ~w", [SG, SH, SM])
    ).
rewrite(h(Q), [M, h(Q)|T], [M2|T], Msg) :-
    hconj(Q, M, M2), !,
    gate_string(M, S1), gate_string(M2, S2),
    format(string(Msg), "H·~w·H -> ~w", [S1, S2]).

step(Gs, Gs2, Msg) :-
    append(Pre, [G|Rest], Gs),
    rewrite(G, Rest, Rest2, Msg), !,
    append(Pre, Rest2, Gs2).

%!  simplify_gates(+Gates, -Simplified, -Log) is det.
simplify_gates(Gs, Out, [M|Ms]) :-
    step(Gs, Gs1, M), !,
    simplify_gates(Gs1, Out, Ms).
simplify_gates(Gs, Gs, []).

%!  approximate(+Gates, +EpsRadians, -Kept, -Dropped)
approximate(Gs, Eps, Kept, Dropped) :-
    partition([G]>>( G = p(_, T), mod1(T, R), abs(2*pi*R) < Eps ), Gs, Dropped, Kept).

% ------------------------------------------------------------- lowering
%!  lower_to_cx(+Gates, -Lowered, -Log) is det.
lower_to_cx([], [], []).
lower_to_cx([p(Q, T)|R0], Out, Log) :- !,
    take_phases([p(Q, T)|R0], Block, Rest),
    lower_block(Block, B2, L1),
    lower_to_cx(Rest, R2, L2),
    append(B2, R2, Out), append(L1, L2, Log).
lower_to_cx([G|T], [G|Out], Log) :- lower_to_cx(T, Out, Log).

take_phases([p(Q, T)|R], [p(Q, T)|B], Rest) :- !, take_phases(R, B, Rest).
take_phases(R, [], R).

lower_block(B, B, []) :- forall(member(p(Qs, _), B), Qs = [_]), !.
lower_block(B, Gates, [L1, L2]) :-
    findall(M-T, ( member(p(Qs, T), B), findall(x(Q), member(Q, Qs), M0), msort(M0, M) ), P0),
    poly_norm(P0, P),
    poly_string(P, PS),
    maplist(gate_string, B, BS), atomic_list_concat(BS, ' ', BA),
    format(string(L1), "diagonal block ~w  =>  φ = ~w", [BA, PS]),
    to_parity(P, Par),
    reduce_mod1(Par, Red, Can, _),
    length(Par, NP), length(Can, NC), length(Red, NR),
    findall(S, ( member(Ts-C, Red), coef_string(C, CS),
                 findall(QQ, member(x(QQ), Ts), Qq), atomic_list_concat(Qq, '⊕', QA),
                 format(string(S), "~w·[~w]", [CS, QA]) ), Ss),
    atomic_list_concat(Ss, ', ', SA),
    format(string(L2), "  parity basis: ~d parities, ~d cancel mod 1, ~d remain: ~w",
           [NP, NC, NR, SA]),
    findall(Key-par(Tgt, Ctrls, C),
            ( member(Ts-C, Red),
              findall(QQ, member(x(QQ), Ts), Qs0), msort(Qs0, Qs1),
              last(Qs1, Tgt), append(Ctrls, [Tgt], Qs1),
              foldl([Qg, Mg0, Mg1]>>(Mg1 is Mg0 \/ (1 << Qg)), Ctrls, 0, Mask),
              gray_rank(Mask, Rank), Key = Tgt-Rank ),
            Keyed),
    keysort(Keyed, Sorted), pairs_values(Sorted, Pars),
    foldl(emit_parity, Pars, [], RevGates),
    reverse(RevGates, Gates).

emit_parity(par(Tgt, Ctrls, C), Acc0, Acc) :-
    findall(cx(Ct, Tgt), member(Ct, Ctrls), Ladder),
    reverse(Ladder, Undo),
    append(Ladder, [p([Tgt], C)|Undo], Seq),
    reverse(Seq, RSeq),
    append(RSeq, Acc0, Acc).

gray_rank(0, 0) :- !.
gray_rank(M, R) :- gray_rank_(M, 0, R).
gray_rank_(0, B, B) :- !.
gray_rank_(M, B0, B) :- B1 is B0 xor M, M1 is M >> 1, gray_rank_(M1, B1, B).
