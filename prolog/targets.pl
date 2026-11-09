:- encoding(utf8).
:- module(targets, [target/3, target_info/4]).
/** <module> Targets from their mathematical definitions

Matrices are built straight from each definition (never from a circuit), so
every derived circuit is checked against first principles.  Matrix = list
of columns; column k lists U[j,k] for j = 0..N-1.
*/
:- use_module(library(lists)).
:- use_module(sim).

%!  target_info(?Name, -DefaultN, -FixedSize, -Description)
target_info(qft,       3, false, "quantum Fourier transform").
target_info(iqft,      3, false, "inverse QFT").
target_info(hadamard,  3, false, "Walsh-Hadamard transform H^n").
target_info(chirp,     3, false, "quadratic phase e^{iπk²/N}").
target_info(diffusion, 3, false, "Grover diffusion 2|s><s| - I").
target_info(mcz,       3, false, "multi-controlled Z").
target_info(cx,        2, true,  "CNOT (control q0, target q1)").
target_info(cz,        2, true,  "controlled-Z").
target_info(swap,      2, true,  "SWAP").
target_info(iswap,     2, true,  "iSWAP").
target_info(ccz,       3, true,  "doubly-controlled Z").
target_info(toffoli,   3, true,  "Toffoli (controls q0,q1, target q2)").
target_info(fredkin,   3, true,  "Fredkin (control q0 swaps q1,q2)").

%!  target(+Name, +N, -Matrix)
target(Name, N, Cols) :-
    Len is 1 << N, Hi is Len - 1,
    findall(Col, ( between(0, Hi, K),
                   findall(Z, ( between(0, Hi, J), el(Name, N, J, K, Z) ), Col) ),
            Cols).

bit(X, Q, B) :- B is (X >> Q) /\ 1.

flat(N, Turns, c(A, B)) :- cexp(Turns, c(C, S)), R is 1/sqrt(1 << N), A is R*C, B is R*S.
diag(J, K, Turns, Z) :- ( J =:= K -> cexp(Turns, Z) ; Z = c(0.0, 0.0) ).
perm(J, Fk, c(V, 0.0)) :- ( J =:= Fk -> V = 1.0 ; V = 0.0 ).

% U[j,k] = e^{2πi jk/N}/√N, exact rational phase
el(qft, N, J, K, Z) :- T is (J*K) rdiv (1 << N), flat(N, T, Z).
el(iqft, N, J, K, Z) :- T is -(J*K) rdiv (1 << N), flat(N, T, Z).
el(hadamard, N, J, K, Z) :- T is (popcount(J /\ K) mod 2) rdiv 2, flat(N, T, Z).
el(chirp, N, J, K, Z) :- T is (K*K) rdiv (2 * (1 << N)), diag(J, K, T, Z).
el(diffusion, N, J, K, c(V, 0.0)) :- D is 2.0/(1 << N), ( J =:= K -> V is D - 1 ; V = D ).
el(mcz, N, J, K, Z) :- Hi is (1 << N) - 1, ( K =:= Hi -> T is 1 rdiv 2 ; T = 0 ), diag(J, K, T, Z).
el(cz, _, J, K, Z) :- ( K =:= 3 -> T is 1 rdiv 2 ; T = 0 ), diag(J, K, T, Z).
el(ccz, _, J, K, Z) :- ( K =:= 7 -> T is 1 rdiv 2 ; T = 0 ), diag(J, K, T, Z).
el(cx, _, J, K, Z) :- bit(K, 0, B), F is K xor (B << 1), perm(J, F, Z).
el(toffoli, _, J, K, Z) :- bit(K, 0, A), bit(K, 1, B), F is K xor ((A /\ B) << 2), perm(J, F, Z).
el(swap, _, J, K, Z) :- bit(K, 0, A), bit(K, 1, B), F is (A << 1) \/ B, perm(J, F, Z).
el(fredkin, _, J, K, Z) :-
    bit(K, 0, C), bit(K, 1, A), bit(K, 2, B),
    (   C =:= 0 -> F = K ; F is 1 \/ (A << 2) \/ (B << 1) ),
    perm(J, F, Z).
el(iswap, _, J, K, Z) :-
    (   J =:= K, ( K =:= 0 ; K =:= 3 ) -> Z = c(1.0, 0.0)
    ;   J + K =:= 3, J * K =:= 2 -> Z = c(0.0, 1.0)
    ;   Z = c(0.0, 0.0)
    ).
