:- encoding(utf8).
:- module(poly, [poly_norm/2, poly_add/3, poly_mul/3, poly_scale/3, poly_const/2,
                 poly_var/2, register/3, mod1/2, is_zero/1, snap/2, reduce_mod1/4,
                 to_parity/2, mono_string/2, coef_string/2, poly_string/2, split_jk/3]).
/** <module> Boolean phase polynomials

A phase polynomial is a multilinear polynomial over Boolean variables
(x*x = x) whose value is a phase in *turns* (1 turn = 2*pi rad).

  Poly  = list of Mono-Coef, sorted by Mono, no zero coefficients
  Mono  = ordset of variables j(A), k(B) or x(Q)
  Coef  = rational (exact, via rdiv) or float

Because exp(2*pi*i*m) = 1 for every integer m, coefficients only matter
mod 1, and every term whose coefficient is an integer CANCELS.  This is the
step that turns the n^2 terms of j*k/N into the gates of the QFT.
*/
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(ordsets)).
:- use_module(library(pairs)).

%!  poly_norm(+Raw, -Poly) is det.
%   Sort, add like terms, drop zeros.
poly_norm(P0, P) :-
    keysort(P0, S),
    sum_runs(S, P).

sum_runs([], []).
sum_runs([M-C|T], Out) :-
    take_same(M, T, C, Sum, Rest),
    (   Sum =:= 0
    ->  Out = Out1
    ;   Out = [M-Sum|Out1]
    ),
    sum_runs(Rest, Out1).

take_same(M, [M2-C2|T], Acc, Sum, Rest) :-
    M2 == M, !,
    Acc1 is Acc + C2,
    take_same(M, T, Acc1, Sum, Rest).
take_same(_, Rest, Sum, Sum, Rest).

poly_add(A, B, C) :-
    append(A, B, AB),
    poly_norm(AB, C).

poly_scale(P, K, Q) :-
    findall(M-C2, (member(M-C, P), C2 is C*K), Q0),
    poly_norm(Q0, Q).

%   Boolean product: monomials multiply by set union (x*x = x).
poly_mul(A, B, C) :-
    findall(M-Cm,
            ( member(M1-C1, A), member(M2-C2, B),
              ord_union(M1, M2, M), Cm is C1*C2 ),
            L),
    poly_norm(L, C).

poly_const(C, P) :- poly_norm([[]-C], P).
poly_var(V, [[V]-1]).

%!  register(+Name, +N, -Poly) is det.
%   The integer x = sum_a 2^a * x_a as a polynomial in its bits Name(a).
register(Name, N, P) :-
    N1 is N - 1,
    findall([V]-W, ( between(0, N1, A), V =.. [Name, A], W is 2^A ), P).

%!  mod1(+C, -R) is det.   R = C mod 1 in (-1/2, 1/2].
mod1(C, R) :-
    R0 is C - floor(C),
    (   R0 > 1 rdiv 2
    ->  R is R0 - 1
    ;   R = R0
    ).

is_zero(C) :-
    (   rational(C)
    ->  C =:= 0
    ;   abs(C) < 1.0e-9
    ).

%!  snap(+X, -R) is det.
%   Floats from numerics are snapped to a nearby dyadic rational so that all
%   further cancellation is exact.
snap(X, R) :-
    (   rational(X)
    ->  R = X
    ;   D is 2^24,
        R0 is round(X*D) rdiv D,
        (   abs(X - R0) < 1.0e-9
        ->  R = R0
        ;   R = X
        )
    ).

%!  reduce_mod1(+Poly, -Reduced, -Cancelled, -Kept) is det.
%   Cancelled: Mono-Coef with integer coefficient.  Kept: Mono-Coef-Reduced.
reduce_mod1([], [], [], []).
reduce_mod1([M-C|T], Red, Can, Kept) :-
    mod1(C, R),
    (   is_zero(R)
    ->  Red = Red1, Can = [M-C|Can1], Kept = Kept1
    ;   Red = [M-R|Red1], Can = Can1, Kept = [M-C-R|Kept1]
    ),
    reduce_mod1(T, Red1, Can1, Kept1).

%!  to_parity(+Poly, -Parities) is det.
%   prod_{i in S} x_i = 2^(1-|S|) * sum_{T subset S, T \= []} (-1)^(|T|-1) * XOR(T)
%   Result: Poly whose "monomials" are the variable sets of parities.
%   The constant term (global phase) is dropped.
to_parity(P, Par) :-
    findall(T-V,
            ( member(M-C, P), M \== [],
              length(M, L),
              subset_of(M, T), T \== [],
              length(T, LT),
              V is C * (1 rdiv 2^(L-1)) * (-1)^(LT-1) ),
            L0),
    poly_norm(L0, Par).

subset_of([], []).
subset_of([X|Xs], [X|Ys]) :- subset_of(Xs, Ys).
subset_of([_|Xs], Ys) :- subset_of(Xs, Ys).

split_jk(M, J, K) :-
    findall(A, member(j(A), M), J),
    findall(B, member(k(B), M), K).

% ---------------------------------------------------------------- printing
var_string(V, S) :- V =.. [F, A], format(string(S), "~w~w", [F, A]).

mono_string([], "1") :- !.
mono_string(M, S) :-
    maplist(var_string, M, Ss),
    atomic_list_concat(Ss, '·', A),
    atom_string(A, S).

coef_string(C, S) :-
    (   integer(C)
    ->  format(string(S), "~d", [C])
    ;   rational(C, N, D)
    ->  format(string(S), "~d/~d", [N, D])
    ;   format(string(S), "~6g", [C])
    ).

poly_string([], "0") :- !.
poly_string(P, S) :-
    findall(T, ( member(M-C, P), coef_string(C, CS),
                 (   M == []
                 ->  T = CS
                 ;   mono_string(M, MS), format(string(T), "~w·~w", [CS, MS])
                 ) ), Ts),
    atomic_list_concat(Ts, ' + ', A),
    atom_string(A, S).
