:- encoding(utf8).
/*  Run:  swipl -g run_tests -t halt tests/test_iq.pl   */
:- use_module('../prolog/iq').
:- use_module(library(plunit)).
:- use_module(library(lists)).

exact(derivation(_, _, _, _, _, E)) :- ( E == proved -> true ; number(E), E < 1.0e-7 ).
proved(derivation(_, _, _, _, _, E)) :- E == proved.
gates(derivation(_, _, _, circ(_, Gs, _, _), _, _), Gs).
labels(D, Ls) :- gates(D, Gs), maplist(gate_label, Gs, Ls0), msort(Ls0, Ls).

:- begin_tests(phase_polynomial).

test(qft_terms_cancel) :-          % jk/16: 16 terms, 6 integer ones cancel
    register(j, 4, J), register(k, 4, K), poly_mul(J, K, JK),
    poly_scale(JK, 1 rdiv 16, P), reduce_mod1(P, R, C, _),
    length(P, 16), length(C, 6), length(R, 10).

test(parity_identity) :-           % x0*x1 = (x0 + x1 - x0^x1)/2
    to_parity([[x(0), x(1)]-1], Par),
    Par = [[x(0)]-A, [x(0), x(1)]-B, [x(1)]-C],
    A =:= 1 rdiv 2, B =:= -1 rdiv 2, C =:= 1 rdiv 2.

:- end_tests(phase_polynomial).

:- begin_tests(derive).

test(qft, [forall(between(1, 5, N))]) :-
    derive(qft, N, [], D), exact(D), gates(D, Gs), length(Gs, L),
    L =:= N + N*(N-1)//2 + N//2.
test(iqft) :- derive(iqft, 4, [], D), exact(D).
test(qft_no_swaps) :- derive(qft, 4, [swaps(false)], D), exact(D),
    labels(D, Ls), \+ memberchk(swap, Ls).
test(qft_lowered_to_cx) :- derive(qft, 3, [basis(cx)], D), exact(D),
    labels(D, Ls), \+ memberchk(cp, Ls).
test(approx_qft) :-
    derive(qft, 6, [approx(0.2)], D),
    D = derivation(_, _, _, C, _, _), target(qft, 6, U), unitary(C, V),
    fidelity(U, V, F), F > 0.95, gates(D, Gs), length(Gs, L), L < 24.
test(hadamard) :- derive(hadamard, 3, [], D), exact(D), labels(D, [h, h, h]).
test(cz_to_cx) :- derive(cz, _, [basis(cx)], D), exact(D),
    labels(D, Ls), msort(Ls, S), S = [cx, cx, s, s, sdg].
test(toffoli) :- derive(toffoli, _, [], D), exact(D), labels(D, [ccz, h, h]).
test(toffoli_clifford_t) :- derive(toffoli, _, [basis(cx)], D), exact(D),
    labels(D, Ls), include(==(cx), Ls, Cx), length(Cx, NC), NC =< 8.
test(cx) :- derive(cx, _, [], D), exact(D), labels(D, [cz, h, h]).
test(swap) :- derive(swap, _, [], D), exact(D), labels(D, [cx, cx, cx]).
test(fredkin) :- derive(fredkin, _, [], D), exact(D), gates(D, Gs), length(Gs, 3).
test(diffusion) :- derive(diffusion, 3, [], D), exact(D).
test(chirp) :- derive(chirp, 3, [], D), exact(D).
test(mcz) :- derive(mcz, 4, [], D), exact(D), gates(D, [_]).
test(iswap) :- derive(iswap, _, [], D), exact(D).
test(symbolic_qft_12) :- derive(qft, 12, [], D), gates(D, Gs), length(Gs, L), L =:= 12 + 66 + 6.

:- end_tests(derive).

:- begin_tests(prove).

test(all_targets_proved, [forall(member(T-N-O, [qft-4-[], iqft-3-[], qft-4-[swaps(false)],
        qft-3-[basis(cx)], hadamard-3-[], cz-2-[basis(cx)], ccz-3-[basis(cx)],
        toffoli-3-[], toffoli-3-[basis(cx)], cx-2-[], swap-2-[], fredkin-3-[],
        iswap-2-[], chirp-3-[], diffusion-3-[], mcz-4-[]]))]) :-
    derive(T, N, O, D), proved(D).
test(large_qft_proved) :- derive(qft, 16, [], D), proved(D).
test(approx_not_proved) :- derive(qft, 6, [approx(0.2)], D), \+ proved(D).
test(qft3_proof_steps) :-
    derive(qft, 3, [], D), D = derivation(_, _, _, C, _, _),
    prove(C, qft, true, Log),
    include([L]>>sub_string(L, _, _, _, "[HH]"), Log, HH), length(HH, 3).

:- end_tests(prove).

:- begin_tests(simplify).

test(cancel_hh) :- simplify_gates([h(0), h(0)], [], _).
test(cancel_through_commuting) :-     % phase on the control commutes with CX
    Q is 1 rdiv 4,
    simplify_gates([cx(0, 1), p([0], Q), cx(0, 1)], G, _), G == [p([0], Q)].
test(merge_phases) :-
    E is 1 rdiv 8,
    simplify_gates([p([0, 1], E), h(2), p([0, 1], E)], G, _),
    G = [h(2), p([0, 1], T)], T =:= 1 rdiv 4.
test(hxh) :- simplify_gates([h(0), x(0), h(0)], [p([0], T)], _), T =:= 1 rdiv 2.
test(no_cancel_across_blocker) :-
    simplify_gates([h(0), cx(0, 1), h(0)], G, _), length(G, 3).

:- end_tests(simplify).
