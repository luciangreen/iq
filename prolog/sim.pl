:- encoding(utf8).
:- module(sim, [apply_gate/4, run_gates/4, basis_state/3, unitary/2, mat_mul/3,
                distance/3, fidelity/3, cexp/2, cmul/3, cadd/3, cabs/2, carg_turns/2,
                identity/2, entry/4]).
/** <module> State-vector simulation

Complex numbers are c(Re, Im) floats.  A state on n qubits is a term
s(A0, ..., A_{N-1}).  A matrix is a list of columns, each a list of c/2.
A circuit is circ(N, Gates, GlobalPhaseTurns, OutputMap) where OutputMap is
`none` or a list giving, for each output bit a, the qubit that holds it.
*/
:- use_module(library(lists)).
:- use_module(library(apply)).

cadd(c(A, B), c(C, D), c(E, F)) :- E is A + C, F is B + D.
csub(c(A, B), c(C, D), c(E, F)) :- E is A - C, F is B - D.
cmul(c(A, B), c(C, D), c(E, F)) :- E is A*C - B*D, F is A*D + B*C.
cconj(c(A, B), c(A, NB)) :- NB is -B.
cscale(K, c(A, B), c(C, D)) :- C is K*A, D is K*B.
cabs(c(A, B), R) :- R is sqrt(A*A + B*B).
cexp(T, c(C, S)) :- X is 2*pi*T, C is cos(X), S is sin(X).
carg_turns(c(A, B), T) :- T is atan2(B, A) / (2*pi).

amp(S, I, A) :- I1 is I + 1, arg(I1, S, A).

%!  apply_gate(+Gate, +N, +S0, -S) is det.
apply_gate(G, N, S0, S) :-
    Len is 1 << N, Hi is Len - 1,
    findall(A, ( between(0, Hi, I), out_amp(G, S0, I, A) ), L),
    S =.. [s|L].

out_amp(h(Q), S, I, A) :-
    M is 1 << Q, I0 is I /\ \ M, I1 is I \/ M,
    amp(S, I0, A0), amp(S, I1, A1), R is 1/sqrt(2),
    (   I /\ M =:= 0
    ->  cadd(A0, A1, B)
    ;   csub(A0, A1, B)
    ),
    cscale(R, B, A).
out_amp(x(Q), S, I, A) :-
    J is I xor (1 << Q), amp(S, J, A).
out_amp(p(Qs, T), S, I, A) :-
    foldl([Q, M0, M1]>>(M1 is M0 \/ (1 << Q)), Qs, 0, Mask),
    amp(S, I, A0),
    (   I /\ Mask =:= Mask
    ->  cexp(T, E), cmul(A0, E, A)
    ;   A = A0
    ).
out_amp(cx(C, T), S, I, A) :-
    (   I /\ (1 << C) =\= 0 -> J is I xor (1 << T) ; J = I ),
    amp(S, J, A).
out_amp(ccx(C1, C2, T), S, I, A) :-
    (   I /\ (1 << C1) =\= 0, I /\ (1 << C2) =\= 0
    ->  J is I xor (1 << T)
    ;   J = I
    ),
    amp(S, J, A).
out_amp(swap(P, Q), S, I, A) :-
    BP is (I >> P) /\ 1, BQ is (I >> Q) /\ 1,
    (   BP =:= BQ -> J = I ; J is I xor ((1 << P) \/ (1 << Q)) ),
    amp(S, J, A).

run_gates(Gs, N, S0, S) :- foldl([G, A, B]>>apply_gate(G, N, A, B), Gs, S0, S).

basis_state(N, X, S) :-
    Len is 1 << N, Hi is Len - 1,
    findall(A, ( between(0, Hi, I), ( I =:= X -> A = c(1.0, 0.0) ; A = c(0.0, 0.0) ) ), L),
    S =.. [s|L].

identity(N, U) :- unitary(circ(N, [], 0, none), U).

%!  unitary(+Circ, -Columns) is det.
unitary(circ(N, Gs, GP, Map), Cols) :-
    Len is 1 << N, Hi is Len - 1,
    cexp(GP, Ph),
    findall(Col,
            ( between(0, Hi, X),
              basis_state(N, X, S0),
              run_gates(Gs, N, S0, S1),
              S1 =.. [s|L1],
              maplist(cmul(Ph), L1, L2),
              relabel(Map, N, L2, Col) ),
            Cols).

relabel(none, _, L, L) :- !.
relabel(Map, N, L, Out) :-
    S =.. [s|L],
    Len is 1 << N, Hi is Len - 1,
    length(Map, NB), NB1 is NB - 1, numlist(0, NB1, Bits),
    findall(A, ( between(0, Hi, Y),
                 foldl([Q, B, X0, X1]>>(X1 is X0 \/ (((Y >> B) /\ 1) << Q)),
                       Map, Bits, 0, X),
                 amp(S, X, A) ),
            Out).

%!  entry(+Cols, +Row, +Col, -Z)
entry(Cols, J, K, Z) :- nth0(K, Cols, C), nth0(J, C, Z).

%!  mat_mul(+A, +B, -C)   C = A*B, all as column lists
mat_mul(A, B, C) :- maplist(mat_vec(A), B, C).
mat_vec(A, V, W) :-
    A = [C0|_], length(C0, L),
    length(Z, L), maplist(=(c(0.0, 0.0)), Z),
    foldl(axpy, A, V, Z, W).
axpy(Col, K, Acc0, Acc) :- maplist([X, Y, R]>>(cmul(K, X, KX), cadd(Y, KX, R)), Col, Acc0, Acc).

flat_cols(Cols, L) :- append(Cols, L).

%!  distance(+U, +V, -D)   ||U - e^{ig} V|| with optimal global phase g
distance(U, V, D) :-
    flat_cols(U, Lu), flat_cols(V, Lv),
    foldl([Xa, Ya, Aa0, Aa1]>>(cconj(Ya, Yac), cmul(Yac, Xa, Pa), cadd(Aa0, Pa, Aa1)), Lu, Lv, c(0.0, 0.0), IP),
    cabs(IP, R),
    (   R > 1.0e-12 -> IP = c(A, B), Ph = c(A/R, B/R) ; Ph = c(1.0, 0.0) ),
    Ph = c(PA, PB), PA1 is PA, PB1 is PB,
    foldl([Xb, Yb, Sb0, Sb1]>>(cmul(c(PA1, PB1), Yb, PYb), csub(Xb, PYb, Dfb), cabs(Dfb, Mb), Sb1 is Sb0 + Mb*Mb),
          Lu, Lv, 0.0, SS),
    D is sqrt(SS).

%!  fidelity(+U, +V, -F)   |Tr(U^dagger V)| / N
fidelity(U, V, F) :-
    flat_cols(U, Lu), flat_cols(V, Lv), length(U, N),
    foldl([Xf, Yf, Af0, Af1]>>(cconj(Xf, Xfc), cmul(Xfc, Yf, Pf), cadd(Af0, Pf, Af1)), Lu, Lv, c(0.0, 0.0), IP),
    cabs(IP, R), F is R / N.
