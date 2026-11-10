:- encoding(utf8).
/*  Checks iq against the 200 input/output pairs in data/io_pairs.pl.

    Run:  swipl -g run_tests -t halt tests/test_io_pairs.pl

    For every (target, n) in the data, iq derives the circuit, simulates each
    listed basis input, and compares with the expected output state up to one
    global phase per circuit.  A second test checks the pairs against the
    definitions in targets.pl directly.
*/
:- use_module('../prolog/iq').
:- use_module(library(plunit)).
:- use_module(library(lists)).
:- use_module(library(apply)).
:- consult('../data/io_pairs.pl').

%   Expected output state as a full list of c(Re, Im).
expected(N, Out, Col) :-
    Len is 1 << N, Hi is Len - 1,
    findall(Z, ( between(0, Hi, J),
                 (   memberchk(J-amp(P, T), Out)
                 ->  R is sqrt(P), cexp(T, c(C, S)), A is R*C, B is R*S, Z = c(A, B)
                 ;   Z = c(0.0, 0.0)
                 ) ),
            Col).

cases(Cases) :- findall(T-N, io(_, T, N, _, _), L), sort(L, Cases).

%   Matrices built from the pairs, column by column, in input order.
pair_matrix(T, N, Ks, Cols) :-
    findall(K-Col, ( io(_, T, N, K, Out), expected(N, Out, Col) ), KC0),
    keysort(KC0, KC), pairs_keys_values(KC, Ks, Cols).

select_cols(Ks, U, Cols) :- findall(C, ( member(K, Ks), nth0(K, U, C) ), Cols).

:- begin_tests(io_pairs).

test(count) :- aggregate_all(count, io(_, _, _, _, _), 200).

test(normalised, [forall(io(Id, _, _, _, Out))]) :-
    foldl([_-amp(P, _), A0, A]>>(A is A0 + P), Out, 0, Sum),
    (   Sum =:= 1 -> true ; format(user_error, "pair ~w not normalised~n", [Id]), fail ).

test(matches_definition, [forall((cases(Cs), member(T-N, Cs)))]) :-
    pair_matrix(T, N, Ks, Cols),
    target(T, N, U), select_cols(Ks, U, Def),
    distance(Cols, Def, D), D < 1.0e-9.

test(matches_derived_circuit, [forall((cases(Cs), member(T-N, Cs)))]) :-
    pair_matrix(T, N, Ks, Cols),
    derive(T, N, [], derivation(_, _, _, Circ, _, _)),
    unitary(Circ, U), select_cols(Ks, U, Got),
    distance(Cols, Got, D),
    (   D < 1.0e-9 -> true
    ;   format(user_error, "~w on ~d qubits: distance ~e~n", [T, N, D]), fail
    ).

:- end_tests(io_pairs).
