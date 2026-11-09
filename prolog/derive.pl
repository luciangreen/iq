:- encoding(utf8).
:- module(derive, [derive/4, derive_qft/5, derive_matrix/5, print_derivation/2]).
/** <module> First-principles derivation

derive(+Name, +N, +Options, -Derivation)

Derivation = derivation(Name, Method, RawCirc, Circ, Trace, Error)

Options:
  basis(native|cx)   cx: lower diagonal blocks to CNOT + phase (parity basis)
  simplify(Bool)     default true
  approx(Eps)        drop phase rotations below Eps radians
  swaps(Bool)        default true; false leaves the output qubits permuted
  verify(Bool)       default: true for n =< 6

Pipeline for a target U:
  1. diagonal in some Hadamard frame H_S U H_S ?  -> phase polynomial over k
     (with a fixed-polarity X frame that minimises the term count)
  2. a permutation (classical reversible)?        -> reversible search
  3. flat |U[j,k]| = 1/sqrt N in some frame?       -> phase polynomial over j,k
  4. otherwise                                    -> breadth-first search
The QFT is also derived symbolically from j*k/N with no matrix at all.
*/
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(option)).
:- use_module(poly).
:- use_module(gates).
:- use_module(sim).
:- use_module(synth).
:- use_module(simplify).
:- use_module(search).
:- use_module(targets).
:- use_module(output).

line(Fmt, Args, S) :- format(string(S), Fmt, Args).

derive(Name, N0, Opts, D) :-
    target_info(Name, DN, Fixed, _), !,
    (   ( Fixed == true ; var(N0) ) -> N = DN ; N = N0 ),
    (   Name == qft  -> derive_qft(N, false, 0, Opts, D)
    ;   Name == iqft -> derive_qft(N, true, 0, Opts, D)
    ;   target(Name, N, U),
        derive_matrix(Name, U, N, Opts, D)
    ).
derive(Name, _, _, _) :-
    throw(error(existence_error(target, Name), _)).

% ------------------------------------------------------------ symbolic QFT
derive_qft(N, Inverse, Shift, Opts, D) :-
    Big is 1 << N,
    (   Inverse == true -> Sg = -1, SS = "-", Nm = iqft ; Sg = 1, SS = "", Nm = qft ),
    line("~wQFT on ~d qubits from its definition:", [SS, N], T0a),
    line("  U|k⟩ = N^-½ Σ_j e^{2πi·~wj(k+~w)/N} |j⟩,  N = ~d", [SS, Shift, Big], T0b),
    line("  phase in turns φ(j,k) = ~wj·(k+~w)/~d;  write j = Σ 2^a j_a, k = Σ 2^b k_b",
         [SS, Shift, Big], T0c),
    register(j, N, J), register(k, N, K0),
    poly_const(Shift, SP), poly_add(K0, SP, K),
    poly_mul(J, K, JK),
    Sc is Sg rdiv Big,
    poly_scale(JK, Sc, Phi),
    reduce_with_trace(Phi, Red, T1),
    option(swaps(Swaps), Opts, true),
    synth_flat(Red, N, Swaps, Raw, T2),
    (   option(verify(V), Opts) -> true ; ( N =< 6 -> V = true ; V = false ) ),
    (   V == true
    ->  ( Inverse == true -> TName = iqft ; TName = qft ),
        (   Shift =:= 0
        ->  target(TName, N, Target)
        ;   Target = none
        )
    ;   Target = none
    ),
    append([[T0a, T0b, T0c], T1, ["Read off the circuit:"], T2], Tr),
    atom_concat(Nm, N, Name),
    finish(Name, "symbolic phase polynomial", Raw, Target, Opts, Tr, D).

% ------------------------------------------------------------ any matrix
subsets_by_size(N, Subsets) :-
    N1 is N - 1, numlist(0, N1, Qs),
    findall(S, ( between(0, N, R), length(S, R), subseq_(Qs, S) ), Subsets).

subseq_([], []).
subseq_([X|Xs], [X|Ys]) :- subseq_(Xs, Ys).
subseq_([_|Xs], Ys) :- subseq_(Xs, Ys).

h_frame(N, S, M) :-
    findall(h(Q), member(Q, S), Gs),
    unitary(circ(N, Gs, 0, none), M).

frames(U, N, Frames) :-
    subsets_by_size(N, Subs),
    findall(S-V, ( member(S, Subs), h_frame(N, S, H),
                   mat_mul(H, U, HU), mat_mul(HU, H, V) ), Frames).

is_diag(V) :-
    forall(( nth0(K, V, Col), nth0(J, Col, Z), J =\= K ), ( cabs(Z, R), R < 1.0e-9 )).

is_flat(V, N) :-
    E is 1/sqrt(1 << N),
    forall(( member(Col, V), member(Z, Col) ), ( cabs(Z, R), abs(R - E) < 1.0e-9 )).

is_perm(V) :-
    forall(( member(Col, V), member(c(A, B), Col) ),
           ( abs(B) < 1.0e-9, ( abs(A) < 1.0e-9 ; abs(A - 1) < 1.0e-9 ) )).

diag_entries(V, D) :- findall(Z, ( nth0(K, V, Col), nth0(K, Col, Z) ), D).

%   Möbius transform: exact multilinear expansion of a table of phases.
moebius(L, Bits, Out) :- moebius_(0, Bits, L, Out).
moebius_(I, Bits, L, L) :- I >= Bits, !.
moebius_(I, Bits, L0, L) :-
    Arr =.. [a|L0], length(L0, Len), Hi is Len - 1, Mask is 1 << I,
    findall(V, ( between(0, Hi, Idx), Ix is Idx + 1, arg(Ix, Arr, X),
                 (   Idx /\ Mask =\= 0
                 ->  Jx is (Idx xor Mask) + 1, arg(Jx, Arr, Y), V is X - Y
                 ;   V = X
                 ) ),
            L1),
    I1 is I + 1,
    moebius_(I1, Bits, L1, L).

table_poly(Table, Bits, VarOf, Poly) :-
    moebius(Table, Bits, Co),
    findall(M-C, ( nth0(Idx, Co, C0), snap(C0, C), \+ is_zero(C),
                   B1 is Bits - 1,
                   findall(Var, ( between(0, B1, B), Idx /\ (1 << B) =\= 0, call(VarOf, B, Var) ), M0),
                   msort(M0, M) ),
            P0),
    poly_norm(P0, Poly).

kvar(B, k(B)).
jkvar(N, B, V) :- ( B < N -> V = k(B) ; A is B - N, V = j(A) ).

poly_of_diagonal(D, N, Poly) :-
    maplist(carg_turns, D, Phi),
    table_poly(Phi, N, kvar, Poly).

%   index = j*N + k : low bits are k, high bits are j
poly_of_transform(V, N, Poly) :-
    Len is 1 << N, Hi is Len - 1, Sc is sqrt(Len),
    findall(T, ( between(0, Hi, J), between(0, Hi, K),
                 entry(V, J, K, c(A, B)), A1 is A*Sc, B1 is B*Sc, carg_turns(c(A1, B1), T) ),
            Phi),
    Bits is 2*N,
    table_poly(Phi, Bits, jkvar(N), Poly).

nonconst_terms(P, NT) :- exclude([M-_]>>(M == []), P, L), length(L, NT).

best_polarity(D, N, Best) :-
    subsets_by_size(N, Subs),
    Len is 1 << N, Hi is Len - 1,
    findall(NT-X-DD,
            ( member(X, Subs),
              foldl([Q, M0, M1]>>(M1 is M0 \/ (1 << Q)), X, 0, Mask),
              findall(Z, ( between(0, Hi, I), J is I xor Mask, nth0(J, D, Z) ), DD),
              poly_of_diagonal(DD, N, P), reduce_mod1(P, R, _, _), nonconst_terms(R, NT) ),
            All),
    All = [First|_],
    foldl([E, B0, B]>>( E = NE-_-_, B0 = NB-_-_, ( NE < NB -> B = E ; B = B0 ) ), All, First, Best).

frame_text([], "") :- !.
frame_text(S, T) :- format(string(T), " in the Hadamard frame on qubits ~w (U = H_S·V·H_S)", [S]).

derive_matrix(Name, U, N, Opts, D) :-
    frames(U, N, Frames),
    (   member(S-V, Frames), is_diag(V)
    ->  derive_diagonal(Name, U, N, S, V, Opts, D)
    ;   is_perm(U), reversible_search(Name, U, N, Opts, D)
    ->  true
    ;   member(S-V, Frames), is_flat(V, N),
        catch(derive_flat(Name, U, N, S, V, Opts, D), synthesis_error(_), fail)
    ->  true
    ;   general_search(Name, U, N, Opts, D)
    ).

derive_diagonal(Name, U, N, S, V, Opts, D) :-
    frame_text(S, FT),
    line("~w: diagonal~w: V = diag(e^{2πi φ(k)})", [Name, FT], L0),
    diag_entries(V, Dg),
    best_polarity(Dg, N, NT-X-DD),
    (   X == [] -> LX = []
    ;   line("  fixed-polarity choice: conjugate by X on qubits ~w (literals 1−k_b) => fewest terms (~d)",
             [X, NT], LX1), LX = [LX1]
    ),
    poly_of_diagonal(DD, N, Poly), poly_string(Poly, PS),
    line("  Möbius transform over ~d Boolean variables: φ = ~w", [N, PS], L1),
    reduce_with_trace(Poly, Red, T1),
    synth_diagonal(Red, N, circ(N, Core, GP, none), T2),
    findall(h(Q), member(Q, S), HS), findall(x(Q), member(Q, X), XS),
    append(HS, XS, Pre), reverse(Pre, Post),
    append([Pre, Core, Post], Gs),
    append([[L0], LX, [L1], T1, T2], Tr),
    ( S == [] -> M = "diagonal phase polynomial" ; M = "diagonal phase polynomial + H frame" ),
    finish(Name, M, circ(N, Gs, GP, none), U, Opts, Tr, D).

derive_flat(Name, U, N, S, V, Opts, D) :-
    frame_text(S, FT),
    Big is 1 << N,
    line("~w: all amplitudes 1/√~d~w: V[j,k] = e^{2πi φ(j,k)}/√N", [Name, Big, FT], L0),
    poly_of_transform(V, N, Poly),
    B2 is 2*N,
    line("  Möbius transform over ~d Boolean variables (j and k bits)", [B2], L1),
    reduce_with_trace(Poly, Red, T1),
    option(swaps(Swaps), Opts, true),
    synth_flat(Red, N, Swaps, circ(N, Core, GP, Map), T2),
    findall(h(Q), member(Q, S), HS),
    append([HS, Core, HS], Gs),
    append([[L0, L1], T1, ["Read off the circuit:"], T2], Tr),
    ( S == [] -> M = "Hadamard-type phase polynomial" ; M = "Hadamard-type phase polynomial + H frame" ),
    finish(Name, M, circ(N, Gs, GP, Map), U, Opts, Tr, D).

reversible_search(Name, U, N, Opts, D) :-
    member(Lib, [[x, cx], [x, cx, ccx]]),
    search(U, N, Lib, 6, Gs), !,
    length(Gs, L),
    line("~w: a permutation of basis states (classical reversible); shortest circuit over ~w: ~d gates",
         [Name, Lib, L], T),
    finish(Name, "reversible search", circ(N, Gs, 0, none), U, Opts, [T], D).

general_search(Name, U, N, Opts, D) :-
    option(depth(Max), Opts, 4),
    Lib = [h, s, sdg, t, tdg, x, cx, ccx],
    (   search(U, N, Lib, Max, Gs)
    ->  true
    ;   throw(synthesis_error("no structure recognised and no short circuit found"))
    ),
    length(Gs, L),
    line("~w: no diagonal or Hadamard-type structure in any H frame; breadth-first search found ~d gates",
         [Name, L], T),
    finish(Name, "search", circ(N, Gs, 0, none), U, Opts, [T], D).

% ------------------------------------------------------------ post-processing
finish(Name, Method, Raw, Target, Opts, Tr0, derivation(Name, Method, Raw, Circ, Tr, Err)) :-
    circuit_summary(Raw, RS),
    line("Raw circuit: ~w", [RS], LRaw),
    Raw = circ(N, G0, GP, Map),
    (   option(basis(cx), Opts)
    ->  lower_to_cx(G0, G1, LL),
        prefix_lines(LL, LL2),
        LowT = ["Lower diagonal blocks to CNOT + phase via the parity basis x·y = ½(x + y − x⊕y):"|LL2]
    ;   G1 = G0, LowT = []
    ),
    (   option(simplify(false), Opts)
    ->  G2 = G1, SimT = []
    ;   simplify_gates(G1, G2, SL),
        length(G1, B), length(G2, A),
        line("Simplify: ~d → ~d gates", [B, A], SH),
        prefix_lines(SL, SL2),
        SimT = [SH|SL2]
    ),
    (   option(approx(Eps), Opts)
    ->  approximate(G2, Eps, G3, Dropped), length(Dropped, ND),
        line("Approximate (drop rotations < ~w rad): ~d gates removed", [Eps, ND], AL),
        ApT = [AL]
    ;   G3 = G2, ApT = []
    ),
    Circ = circ(N, G3, GP, Map),
    (   Target == none
    ->  Err = none, VT = []
    ;   unitary(Circ, UC),
        distance(UC, Target, Err),
        fidelity(Target, UC, F),
        ( Err < 1.0e-7 -> Kind = exact ; Kind = approximate ),
        line("Verify against definition: ||U_circ - e^{iγ}U|| = ~e, fidelity ~6f (~w)", [Err, F, Kind], VL),
        VT = [VL]
    ),
    circuit_summary(Circ, FS),
    line("Final circuit: ~w", [FS], FL),
    line("Method: ~w", [Method], ML),
    append([Tr0, [LRaw], LowT, SimT, ApT, VT, [FL, ML]], Tr).

prefix_lines(L, L2) :- maplist([S, T]>>string_concat("  ", S, T), L, L2).

print_derivation(derivation(_, _, _, Circ, Tr, _), Opts) :-
    (   option(quiet(true), Opts) -> true ; forall(member(L, Tr), format("~w~n", [L])), nl ),
    circuit_draw(Circ, Pic), format("~w~n", [Pic]),
    (   option(qasm(true), Opts) -> circuit_qasm(Circ, Q), format("~n~w", [Q]) ; true ).
