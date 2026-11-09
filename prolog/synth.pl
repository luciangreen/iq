:- encoding(utf8).
:- module(synth, [synth_flat/5, synth_diagonal/4, reduce_with_trace/3]).
/** <module> Reading a circuit off a reduced phase polynomial

Hadamard-type transform  U|k> = N^-1/2 sum_j e^{2 pi i phi(j,k)} |j>:

  1. coefficient 1/2 on j_a*k_b gives (-1)^{j_a k_b}: a Hadamard on the
     qubit holding k_b turns it into j_a.  These pairs a -> pi(a) must form
     a permutation.
  2. every other term with output bits J and input bits K is a
     (multi-)controlled phase on qubits pi(J) + K, applied after the
     Hadamards on pi(J) and before the Hadamards on K.  These precedence
     constraints are topologically sorted to give the gate order.
  3. output bit j_a ends on qubit pi(a): SWAPs (or a relabelling) fix it.

Diagonal transform  U|k> = e^{2 pi i phi(k)} |k>: each term c*k_a*k_b*...
is directly a phase gate p([a,b,...], c).
*/
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(pairs)).
:- use_module(poly).
:- use_module(gates).

line(Fmt, Args, S) :- format(string(S), Fmt, Args).

%!  reduce_with_trace(+Poly, -Reduced, -Trace) is det.
reduce_with_trace(P, Red, [L1, L2, L3]) :-
    length(P, NT),
    line("  φ has ~d terms; reduce mod 1 (e^{2πi·m} = 1 for integer m):", [NT], L1),
    reduce_mod1(P, Red, Can, _),
    exclude([Mx-_]>>(Mx == []), Can, CanNC),
    length(CanNC, NC),
    (   NC =:= 0
    ->  L2 = "    nothing cancels"
    ;   findall(S, ( member(M-C, CanNC), mono_string(M, MS), coef_string(C, CS),
                     format(string(S), "~w:~w", [MS, CS]) ), Ss),
        (   length(Ss, LS), LS > 24 -> length(Sh, 24), append(Sh, _, Ss), More = ", …" ; Sh = Ss, More = "" ),
        atomic_list_concat(Sh, ', ', SA),
        line("    cancelled ~d integer terms: ~w~w", [NC, SA, More], L2)
    ),
    exclude([My-_]>>(My == []), Red, RedNC), length(RedNC, NR),
    poly_string(Red, RS),
    line("    surviving ~d terms: φ ≡ ~w", [NR, RS], L3).

%!  synth_diagonal(+Poly, +N, -Circ, -Trace) is det.
synth_diagonal(P, N, circ(N, Gates, GP, none), [L]) :-
    (   memberchk([]-GP, P) -> true ; GP = 0 ),
    findall(G, ( member(M-C, P), M \== [],
                 findall(B, member(k(B), M), Qs), mkphase(Qs, C, G) ), Gates),
    length(Gates, NG),
    line("  each surviving term c·k_a·k_b… is a phase e^{2πi c} on qubits a,b,… => ~d phase gates",
         [NG], L).

is_hkernel([j(_), k(_)]-C) :- C1 is C - 1 rdiv 2, mod1(C1, R), is_zero(R).

%!  synth_flat(+Poly, +N, +Swaps, -Circ, -Trace) is det.
%   Throws synthesis_error(Reason) if the polynomial has no such circuit.
synth_flat(P, N, Swaps, circ(N, Gates, GP, Map), Trace) :-
    (   select([]-GP, P, P1) -> true ; GP = 0, P1 = P ),
    partition(is_hkernel, P1, Kern, Others),
    findall(A-B, member([j(A), k(B)]-_, Kern), Pi0),
    msort(Pi0, Pi),
    pairs_keys_values(Pi, As, Bs),
    N1 is N - 1, numlist(0, N1, All),
    (   As == All, msort(Bs, All)
    ->  true
    ;   throw(synthesis_error("Hadamard kernels do not form a permutation"))
    ),
    findall(S, ( member(A-B, Pi), format(string(S), "j~d←H(q~d)", [A, B]) ), PS),
    atomic_list_concat(PS, ', ', PA),
    line("  coefficient ½ on j_a·k_b  =>  (-1)^{j_a k_b}: a Hadamard on the qubit holding k_b turns it into j_a:", [], T1),
    line("    ~w", [PA], T2),
    % precedence edges  pi(a) -> b
    findall(QA-B,
            ( member(M-_, Others), split_jk(M, J, K),
              member(A, J), member(B, K), memberchk(A-QA, Pi),
              (   QA == B
              ->  mono_string(M, MS),
                  format(string(E), "term ~w needs one qubit to hold j and k at once", [MS]),
                  throw(synthesis_error(E))
              ;   true
              ) ),
            Edges0),
    sort(Edges0, Edges),
    topo(All, Edges, Order),
    findall(S, ( member(Q, Order), format(string(S), "H(q~d)", [Q]) ), OS),
    atomic_list_concat(OS, ' → ', OA),
    line("  a term j_a·k_b must act after H(q_π(a)) and before H(q_b) => order ~w", [OA], T3),
    % schedule
    partition([Mz-_]>>(split_jk(Mz, [], _)), Others, Before, After),
    maplist(term_gate(Pi), Before, BG),
    foldl(emit_after(Pi, Order, After), Order, [], RevMid),
    reverse(RevMid, MidGates),
    append(BG, MidGates, Core),
    length(Others, NO), length(Core, NCore), NCP is NCore - N,
    line("  remaining ~d terms become controlled phases => ~d H + ~d phase gates", [NO, N, NCP], T4),
    pairs_values(Pi, Where),
    (   Where == All
    ->  SwapGates = [], Map = none,
        T5 = "  output bits already in place: no SWAPs"
    ;   Swaps == true
    ->  swap_network(0, N, Where, SwapGates), Map = none,
        length(SwapGates, NS),
        line("  output j_a sits on qubit π(a) = ~w => ~d SWAP(s)", [Where, NS], T5)
    ;   SwapGates = [], Map = Where,
        line("  output bits left permuted (no SWAPs): j_a on qubit ~w", [Where], T5)
    ),
    append(Core, SwapGates, Gates),
    Trace = [T1, T2, T3, T4, T5].

term_gate(Pi, M-C, G) :-
    split_jk(M, J, K),
    findall(Q, ( member(A, J), memberchk(A-Q, Pi) ), QJ),
    append(QJ, K, Qs),
    mkphase(Qs, C, G).

%   After H on qubit Q, emit every term whose last-processed output qubit is Q.
emit_after(Pi, Order, Terms, Q, Acc0, Acc) :-
    findall(G, ( member(M-C, Terms), last_h(Pi, Order, M, Q), term_gate(Pi, M-C, G) ), Gs),
    reverse([h(Q)|Gs], R),
    append(R, Acc0, Acc).

last_h(Pi, Order, M, Q) :-
    split_jk(M, J, _),
    findall(Pos-QA, ( member(A, J), memberchk(A-QA, Pi), nth0(Pos, Order, QA) ), L),
    max_member(_-Q0, L),
    Q0 == Q.

%   Kahn's algorithm, highest-numbered ready qubit first (deterministic).
topo([], _, []) :- !.
topo(Rem, Edges, [Q|Order]) :-
    include([X]>>( \+ ( member(P-X, Edges), memberchk(P, Rem) ) ), Rem, Ready),
    (   Ready == []
    ->  throw(synthesis_error("cyclic ordering constraints between controlled phases"))
    ;   max_list(Ready, Q),
        selectchk(Q, Rem, Rem1),
        topo(Rem1, Edges, Order)
    ).

%   Where[a] = qubit currently holding output bit a.
swap_network(A, N, _, []) :- A >= N, !.
swap_network(A, N, Where, Gs) :-
    nth0(A, Where, QA),
    A1 is A + 1,
    (   QA =:= A
    ->  swap_network(A1, N, Where, Gs)
    ;   nth0(B, Where, A),           % bit B currently sits on qubit A
        set_nth0(A, Where, A, W1),
        set_nth0(B, W1, QA, W2),
        Lo is min(A, QA), Hi is max(A, QA),
        Gs = [swap(Lo, Hi)|Gs1],
        swap_network(A1, N, W2, Gs1)
    ).

set_nth0(I, L, V, L2) :- nth0(I, L, _, R), nth0(I, L2, V, R).
