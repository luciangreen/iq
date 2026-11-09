:- encoding(utf8).
:- module(search, [search/5, instantiate/3, ukey/2]).
/** <module> Breadth-first search for a shortest gate sequence

Used when no algebraic structure is recognised and for classical
reversible (permutation) targets.  States are whole unitaries, deduplicated
up to global phase, so the first hit is a shortest circuit.
*/
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(assoc)).
:- use_module(sim).

qubit(N, Q) :- N1 is N - 1, between(0, N1, Q).

%!  instantiate(+Library, +N, -Gates)
%   Library: list of h, x, s, sdg, t, tdg, z, cx, cz, swap, ccx, cp(Turns).
instantiate(Lib, N, Gs) :-
    findall(G, ( member(L, Lib), inst(L, N, G) ), Gs).

inst(h, N, h(Q)) :- qubit(N, Q).
inst(x, N, x(Q)) :- qubit(N, Q).
inst(z, N, p([Q], T)) :- qubit(N, Q), T is 1 rdiv 2.
inst(s, N, p([Q], T)) :- qubit(N, Q), T is 1 rdiv 4.
inst(sdg, N, p([Q], T)) :- qubit(N, Q), T is -1 rdiv 4.
inst(t, N, p([Q], T)) :- qubit(N, Q), T is 1 rdiv 8.
inst(tdg, N, p([Q], T)) :- qubit(N, Q), T is -1 rdiv 8.
inst(cx, N, cx(C, T)) :- qubit(N, C), qubit(N, T), C \== T.
inst(cz, N, p([A, B], T)) :- qubit(N, A), qubit(N, B), A < B, T is 1 rdiv 2.
inst(cp(Tu), N, p([A, B], Tu)) :- qubit(N, A), qubit(N, B), A < B.
inst(swap, N, swap(A, B)) :- qubit(N, A), qubit(N, B), A < B.
inst(ccx, N, ccx(A, B, T)) :- qubit(N, A), qubit(N, B), A < B, qubit(N, T), T \== A, T \== B.

%!  ukey(+U, -Key)   fingerprint up to global phase
ukey(U, Key) :-
    append(U, L),
    first_big(L, Z), cabs(Z, R), Z = c(A, B), cconj_unit(A, B, R, Ph),
    maplist([X, K]>>(cmul(Ph, X, c(P, Q)), RP is round(P*1.0e6), RQ is round(Q*1.0e6), K = RP/RQ), L, Key).

first_big([Z|T], W) :- ( cabs(Z, R), R > 1.0e-6 -> W = Z ; first_big(T, W) ).
cconj_unit(A, B, R, c(X, Y)) :- X is A/R, Y is -B/R.

apply_u(G, N, U, U2) :-
    maplist([Col, Col2]>>( S =.. [s|Col], apply_gate(G, N, S, S2), S2 =.. [s|Col2] ), U, U2).

%!  search(+Target, +N, +Library, +MaxDepth, -Gates) is semidet.
search(Target, N, Lib, Max, Gates) :-
    instantiate(Lib, N, Gs),
    identity(N, I),
    ukey(Target, TK), ukey(I, IK),
    (   IK == TK
    ->  Gates = []
    ;   empty_assoc(V0), put_assoc(IK, V0, t, V1),
        bfs([I-[]], Gs, N, TK, V1, 1, Max, Rev),
        reverse(Rev, Gates)
    ).

bfs(Frontier, Gs, N, TK, V, D, Max, Path) :-
    D =< Max,
    findall(U2-[G|P], ( member(U-P, Frontier), member(G, Gs), apply_u(G, N, U, U2) ), Kids),
    filter(Kids, TK, V, V2, Next, Found),
    (   Found = found(Path0)
    ->  Path = Path0
    ;   Next \== [],
        D1 is D + 1,
        bfs(Next, Gs, N, TK, V2, D1, Max, Path)
    ).

filter([], _, V, V, [], none).
filter([U-P|T], TK, V0, V, Next, Found) :-
    ukey(U, K),
    (   K == TK
    ->  Found = found(P), V = V0, Next = []
    ;   get_assoc(K, V0, _)
    ->  filter(T, TK, V0, V, Next, Found)
    ;   put_assoc(K, V0, t, V1),
        Next = [U-P|Next1],
        filter(T, TK, V1, V, Next1, Found)
    ).
