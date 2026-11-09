:- encoding(utf8).
:- module(output, [circuit_summary/2, circuit_draw/2, circuit_qasm/2, gate_counts/2, depth/2]).
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(pairs)).
:- use_module(gates).

gate_counts(Gs, Counts) :-
    maplist(gate_label, Gs, Ls),
    msort(Ls, S),
    clumped(S, Counts).

depth(circ(N, Gs, _, _), D) :-
    length(L0, N), maplist(=(0), L0),
    foldl(depth_step, Gs, L0, L),
    max_list([0|L], D).

depth_step(G, L0, L) :-
    gate_qubits(G, Qs),
    findall(V, ( member(Q, Qs), nth0(Q, L0, V) ), Vs),
    max_list(Vs, M), M1 is M + 1,
    foldl([Qd, Ad0, Ad]>>(nth0(Qd, Ad0, _, Rd), nth0(Qd, Ad, M1, Rd)), Qs, L0, L).

circuit_summary(circ(N, Gs, GP, Map), S) :-
    length(Gs, NG),
    include([G]>>(gate_qubits(G, Q), length(Q, LQ), LQ >= 2), Gs, Multi),
    length(Multi, NM),
    depth(circ(N, Gs, GP, Map), D),
    gate_counts(Gs, C),
    findall(X, ( member(K-V, C), format(string(X), "~w:~d", [K, V]) ), Xs),
    atomic_list_concat(Xs, ', ', XA),
    format(string(S), "~d gates (2q+: ~d), depth ~d [~w]", [NG, NM, D, XA]).

%!  circuit_qasm(+Circ, -Text)   OpenQASM 3
circuit_qasm(circ(N, Gs, GP, Map), Text) :-
    format(string(H), "OPENQASM 3.0;~ninclude \"stdgates.inc\";~nqubit[~d] q;", [N]),
    (   GP =:= 0 -> GL = [] ; R is 2*pi*GP, format(string(GS), "gphase(~15g);", [R]), GL = [GS] ),
    maplist(qasm_line, Gs, Ls),
    (   Map == none -> ML = [] ; format(string(MS), "// output bit a is on qubit ~w", [Map]), ML = [MS] ),
    append([[H], GL, Ls, ML], All),
    atomic_list_concat(All, '\n', A),
    format(string(Text), "~w~n", [A]).

% ------------------------------------------------------------ ASCII drawing
cells(cx(C, T), [C-"*", T-"X"]) :- !.
cells(ccx(A, B, T), [A-"*", B-"*", T-"X"]) :- !.
cells(swap(A, B), [A-"x", B-"x"]) :- !.
cells(p([Q], T), [Q-L]) :- !, label1(p([Q], T), L).
cells(p(Qs, T), Cs) :- !,
    last(Qs, Last),
    gate_label(p(Qs, T), Lb),
    (   memberchk(Lb, [cz, ccz]) -> L = "Z" ; angle_string(T, AS), format(string(L), "P(~w)", [AS]) ),
    findall(Q-C, ( member(Q, Qs), ( Q == Last -> C = L ; C = "*" ) ), Cs).
cells(G, [Q-L]) :- gate_qubits(G, [Q]), gate_label(G, Lb), upcase_atom(Lb, U), atom_string(U, L).

label1(G, L) :-
    gate_label(G, Lb),
    (   Lb == p
    ->  G = p(_, T), angle_string(T, AS), format(string(L), "P(~w)", [AS])
    ;   upcase_atom(Lb, U), atom_string(U, L)
    ).

circuit_draw(circ(N, Gs, _, _), Text) :-
    N1 is N - 1, numlist(0, N1, Qs),
    findall(R0, ( member(Q, Qs), format(string(R0), "q~d: ", [Q]) ), Heads),
    foldl(draw_col(Qs), Gs, Heads, Rows),
    atomic_list_concat(Rows, '\n', A),
    atom_string(A, Text).

draw_col(Qs, G, Rows0, Rows) :-
    cells(G, Cs),
    findall(Len, ( member(_-C, Cs), string_length(C, Len) ), Lens),
    max_list(Lens, M), W is M + 2,
    gate_qubits(G, GQ), min_list(GQ, Lo), max_list(GQ, Hi),
    maplist(draw_cell(Cs, W, Lo, Hi), Qs, Rows0, Rows).

draw_cell(Cs, W, Lo, Hi, Q, R0, R) :-
    (   memberchk(Q-C, Cs) -> true
    ;   Q > Lo, Q < Hi -> C = "+"
    ;   C = "-"
    ),
    string_length(C, L),
    Left is (W - L) // 2, Right is W - L - Left,
    dashes(Left, DL), dashes(Right, DR),
    string_concat(R0, DL, R1), string_concat(R1, C, R2), string_concat(R2, DR, R).

dashes(K, S) :- length(L, K), maplist(=(0'-), L), string_codes(S, L).
