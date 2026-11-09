:- encoding(utf8).
:- module(prove, [prove/4, circuit_pathsum/2, pathsum_string/2, reduce_pathsum/3,
                  lift/2, gxor/3, gand/3]).
/** <module> Exact equivalence proofs by sum-over-paths

Every circuit over H, X, CX, CCX, SWAP and phase gates acts as a path sum

    U|x> = 2^(-h/2) * sum_{y in {0,1}^h} e^{2 pi i phi(x,y)} |f(x,y)>

with phi a Boolean phase polynomial (exact rational turns, mod 1) and each
output f_q a GF(2) polynomial.  To prove circuit C = U we build the path sum
of C, append the path sum of U^dagger written from U's *definition*, and
rewrite with two sound rules:

  [Elim]  a path variable y occurring nowhere:  sum_y 1 = 2          h -= 2
  [HH]    y not in the outputs, and every phase term containing y has
          coefficient 1/2, so phi = phi' + 1/2*y*R with R = v xor Q for a
          path variable v not in Q:
              sum_y (-1)^(y*R) = 2*[v = Q]
          => drop y, substitute v := Q everywhere                    h -= 2

If the result is the identity (no path variables, f_q = x_q, constant phase)
then C = U exactly, up to the reported global phase.  Each step is an
equation, so the log is a proof.  For the QFT on n qubits the proof is n
[HH] steps and runs in polynomial time.

Representations
  GF(2) poly : ordset of monomials, monomial = ordset of x(Q) / y(I); [] is 1
  phase      : poly.pl polynomial (Mono-Coef)
  path sum   : ps(N, Outs, Phase, H, PathVars, NextIndex)
*/
:- use_module(library(lists)).
:- use_module(library(apply)).
:- use_module(library(ordsets)).
:- use_module(poly).

% ---------------------------------------------------------------- GF(2)
gxor(A, B, C) :- ord_symdiff(A, B, C).

gand(A, B, C) :-
    findall(M, ( member(M1, A), member(M2, B), ord_union(M1, M2, M) ), L),
    msort(L, S),
    cancel_pairs(S, C).

cancel_pairs([], []).
cancel_pairs([X, Y|T], R) :- X == Y, !, cancel_pairs(T, R).
cancel_pairs([X|T], [X|R]) :- cancel_pairs(T, R).

gsubst(V, Q, P, R) :- foldl(gsub_m(V, Q), P, [], R).
gsub_m(V, Q, M, Acc, R) :-
    (   ord_memberchk(V, M)
    ->  ord_del_element(M, V, M1),
        gand([M1], Q, T),
        gxor(Acc, T, R)
    ;   gxor(Acc, [M], R)
    ).

%!  lift(+GF2, -Poly)  arithmetic value: a xor b = a + b - 2ab
lift(G, P) :- foldl(lift_step, G, [], P).
lift_step(M, Acc, P) :-
    MP = [M-1],
    poly_mul(Acc, MP, AM), poly_scale(AM, -2, AM2),
    poly_add(Acc, MP, S), poly_add(S, AM2, P).

%   phase substitution v := Q  (Q a GF(2) polynomial)
psubst(V, Q, Ph, R) :-
    lift(Q, LQ),
    foldl(psub_t(V, LQ), Ph, [], R0),
    poly_norm(R0, R).
psub_t(V, LQ, M-C, Acc, R) :-
    (   ord_memberchk(V, M)
    ->  ord_del_element(M, V, M1),
        poly_mul([M1-C], LQ, T),
        append(Acc, T, R)
    ;   append(Acc, [M-C], R)
    ).

gstring([], "0") :- !.
gstring(G, S) :-
    maplist(mono_string, G, Ss),
    atomic_list_concat(Ss, ' ⊕ ', A),
    atom_string(A, S).

% ---------------------------------------------------------------- path sums
ps_init(N, ps(N, Outs, [], 0, [], 0)) :-
    N1 is N - 1,
    findall([[x(Q)]], between(0, N1, Q), Outs).

add_phase(P, ps(N, O, Ph0, H, Pa, Nx), ps(N, O, Ph, H, Pa, Nx)) :-
    poly_add(Ph0, P, Ph1),
    reduce_mod1(Ph1, Ph, _, _).

fresh(ps(N, O, Ph, H, Pa, Nx), y(Nx), ps(N, O, Ph, H1, Pa1, Nx1)) :-
    H1 is H + 1, Nx1 is Nx + 1,
    append(Pa, [y(Nx)], Pa1).

out(ps(_, O, _, _, _, _), Q, G) :- nth0(Q, O, G).

set_out(Q, G, ps(N, O, Ph, H, Pa, Nx), ps(N, O1, Ph, H, Pa, Nx)) :-
    nth0(Q, O, _, R), nth0(Q, O1, G, R).

phase_of_product(T, Gs, P) :-
    snap(T, T1),
    foldl([G, A0, A]>>(lift(G, L), poly_mul(A0, L, A)), Gs, [[]-T1], P).

apply_ps(h(Q), PS0, PS) :-
    out(PS0, Q, G), lift(G, L),
    fresh(PS0, Y, PS1),
    poly_mul(L, [[Y]-1], LY),
    poly_scale(LY, 1 rdiv 2, P),
    add_phase(P, PS1, PS2),
    set_out(Q, [[Y]], PS2, PS).
apply_ps(x(Q), PS0, PS) :-
    out(PS0, Q, G), gxor(G, [[]], G2), set_out(Q, G2, PS0, PS).
apply_ps(cx(C, T), PS0, PS) :-
    out(PS0, C, GC), out(PS0, T, GT), gxor(GT, GC, G2), set_out(T, G2, PS0, PS).
apply_ps(ccx(A, B, T), PS0, PS) :-
    out(PS0, A, GA), out(PS0, B, GB), out(PS0, T, GT),
    gand(GA, GB, AB), gxor(GT, AB, G2), set_out(T, G2, PS0, PS).
apply_ps(swap(A, B), PS0, PS) :-
    out(PS0, A, GA), out(PS0, B, GB),
    set_out(A, GB, PS0, PS1), set_out(B, GA, PS1, PS).
apply_ps(p(Qs, T), PS0, PS) :-
    findall(G, ( member(Q, Qs), out(PS0, Q, G) ), Gs),
    phase_of_product(T, Gs, P),
    add_phase(P, PS0, PS).

%!  circuit_pathsum(+Circ, -PathSum)
circuit_pathsum(circ(N, Gs, GP, Map), PS) :-
    ps_init(N, PS0),
    snap(GP, GP1), poly_const(GP1, GPP), add_phase(GPP, PS0, PS1),
    foldl(apply_ps, Gs, PS1, PS2),
    (   Map == none
    ->  PS = PS2
    ;   PS2 = ps(N1, O, Ph, H, Pa, Nx),
        findall(G, ( member(Q, Map), nth0(Q, O, G) ), O2),
        PS = ps(N1, O2, Ph, H, Pa, Nx)
    ).

pathsum_string(ps(_, O, Ph, H, Pa, _), S) :-
    maplist([V, VS]>>(mono_string([V], VS)), Pa, PaS),
    atomic_list_concat(PaS, ',', PA),
    poly_string(Ph, PhS),
    findall(OS, ( nth0(Q, O, G), gstring(G, GS), format(string(OS), "q~d:~w", [Q, GS]) ), OSs),
    atomic_list_concat(OSs, ', ', OA),
    format(string(S0), "2^(-~d/2) Σ_{~w} e^{2πi[~w]} |~w⟩", [H, PA, PhS, OA]),
    truncate(S0, 400, S).

truncate(S0, Max, S) :-
    string_length(S0, L),
    (   L =< Max -> S = S0
    ;   sub_string(S0, 0, Max, _, S1), string_concat(S1, " …", S)
    ).

% ---------------------------------------------------------------- reduction
has_var(Y, M-_) :- ord_memberchk(Y, M).
in_outputs(Y, O) :- member(G, O), member(M, G), ord_memberchk(Y, M), !.
half(C) :- C1 is C - 1 rdiv 2, mod1(C1, R), is_zero(R).
rem_y(Y, M-_, Acc, R) :- ord_del_element(M, Y, M1), gxor(Acc, [M1], R).

%!  reduce_pathsum(+PS0, -PS, -Log) is det.
reduce_pathsum(PS0, PS, [M|Ms]) :-
    rule(PS0, PS1, M), !,
    reduce_pathsum(PS1, PS, Ms).
reduce_pathsum(PS, PS, []).

rule(ps(N, O, Ph, H, Pa, Nx), ps(N, O, Ph, H2, Pa2, Nx), Msg) :-          % [Elim]
    member(Y, Pa), \+ in_outputs(Y, O),
    \+ ( member(T, Ph), has_var(Y, T) ), !,
    selectchk(Y, Pa, Pa2), H2 is H - 2,
    mono_string([Y], YS),
    format(string(Msg), "[Elim] ~w occurs nowhere: Σ_~w 1 = 2", [YS, YS]).
rule(ps(N, O, Ph, H, Pa, Nx), ps(N, O2, Ph3, H2, Pa2, Nx), Msg) :-         % [HH]
    member(Y, Pa), \+ in_outputs(Y, O),
    include(has_var(Y), Ph, Ty), Ty \== [],
    forall(member(_-C, Ty), half(C)),
    foldl(rem_y(Y), Ty, [], R),
    member(V, Pa), V \== Y,
    ord_memberchk([V], R),
    \+ ( member(M2, R), M2 \== [V], ord_memberchk(V, M2) ), !,
    gxor(R, [[V]], Q),
    exclude(has_var(Y), Ph, Ph1),
    psubst(V, Q, Ph1, Ph2),
    reduce_mod1(Ph2, Ph3, _, _),
    maplist(gsubst(V, Q), O, O2),
    subtract(Pa, [Y, V], Pa2),
    H2 is H - 2,
    mono_string([Y], YS), mono_string([V], VS), gstring(R, RS), gstring(Q, QS),
    format(string(Msg), "[HH] Σ_~w (-1)^{~w·(~w)} = 2·[~w = ~w]  =>  ~w := ~w",
           [YS, YS, RS, VS, QS, VS, QS]).

is_identity(ps(N, O, Ph, H, Pa, _)) :-
    Pa == [], H =:= 0,
    forall(member(M-_, Ph), M == []),
    N1 is N - 1,
    findall([[x(Q)]], between(0, N1, Q), O).

% ---------------------------------------------------------------- definitions
%   U^dagger for each target, written from its mathematical definition.
register_lift(O, K) :-
    foldl([G, B0-A0, B1-A1]>>(lift(G, L), W is 2^B0, poly_scale(L, W, LW),
                              poly_add(A0, LW, A1), B1 is B0 + 1),
          O, 0-[], _-K).

fresh_n(0, PS, PS, []) :- !.
fresh_n(K, PS0, PS, [Y|Ys]) :- fresh(PS0, Y, PS1), K1 is K - 1, fresh_n(K1, PS1, PS, Ys).

set_outs(O, ps(N, _, Ph, H, Pa, Nx), ps(N, O, Ph, H, Pa, Nx)).

def_inverse(Name, PS0, PS, Desc) :-
    memberchk(Name, [qft, iqft]), !,
    PS0 = ps(N, O, _, _, _, _), Big is 1 << N,
    register_lift(O, K),
    fresh_n(N, PS0, PS1, Zs),
    foldl([Zq, Bq0-Aq0, Bq1-Aq1]>>(Wq is 2^Bq0, poly_add(Aq0, [[Zq]-Wq], Aq1), Bq1 is Bq0 + 1), Zs, 0-[], _-J),
    poly_mul(J, K, JK),
    (   Name == qft -> Sg = -1, SS = "-" ; Sg = 1, SS = "+" ),
    Sc is Sg rdiv Big, poly_scale(JK, Sc, P),
    add_phase(P, PS1, PS2),
    findall([[Z]], member(Z, Zs), Outs), set_outs(Outs, PS2, PS),
    format(string(Desc), "~w†|k⟩ = 2^(-n/2) Σ_z e^{~w2πi·z·k/~d} |z⟩", [Name, SS, Big]).
def_inverse(hadamard, PS0, PS, "H^n†|k⟩ = 2^(-n/2) Σ_z (-1)^{z·k} |z⟩") :- !,
    PS0 = ps(N, O, _, _, _, _),
    fresh_n(N, PS0, PS1, Zs),
    foldl([Zh, Gh, Ah0, Ah1]>>(lift(Gh, Lh), poly_mul([[Zh]-(-1 rdiv 2)], Lh, Th), poly_add(Ah0, Th, Ah1)),
          Zs, O, [], P),
    add_phase(P, PS1, PS2),
    findall([[Z]], member(Z, Zs), Outs), set_outs(Outs, PS2, PS).
def_inverse(Name, PS0, PS, "phase −½·Π_b k_b") :-
    memberchk(Name, [cz, ccz, mcz]), !,
    PS0 = ps(_, O, _, _, _, _),
    phase_of_product(-1 rdiv 2, O, P),
    add_phase(P, PS0, PS).
def_inverse(chirp, PS0, PS, Desc) :- !,
    PS0 = ps(N, O, _, _, _, _), D is 2 * (1 << N),
    register_lift(O, K), poly_mul(K, K, KK),
    poly_scale(KK, -1 rdiv D, P), add_phase(P, PS0, PS),
    format(string(Desc), "phase −k²/~d", [D]).
def_inverse(cx, PS0, PS, "|a,b⟩ → |a, b⊕a⟩") :- !,
    apply_ps(cx(0, 1), PS0, PS).
def_inverse(toffoli, PS0, PS, "|a,b,c⟩ → |a, b, c⊕ab⟩") :- !,
    apply_ps(ccx(0, 1, 2), PS0, PS).
def_inverse(swap, PS0, PS, "|a,b⟩ → |b,a⟩") :- !,
    apply_ps(swap(0, 1), PS0, PS).
def_inverse(fredkin, PS0, PS, "|c,a,b⟩ → |c, a⊕c(a⊕b), b⊕c(a⊕b)⟩") :- !,
    out(PS0, 0, C), out(PS0, 1, A), out(PS0, 2, B),
    gxor(A, B, AB), gand(C, AB, D), gxor(A, D, A2), gxor(B, D, B2),
    set_out(1, A2, PS0, PS1), set_out(2, B2, PS1, PS).
def_inverse(iswap, PS0, PS, "|a,b⟩ → i^{-(a⊕b)} |b,a⟩") :- !,
    out(PS0, 0, A), out(PS0, 1, B), gxor(A, B, AB), lift(AB, L),
    poly_scale(L, -1 rdiv 4, P), add_phase(P, PS0, PS1),
    apply_ps(swap(0, 1), PS1, PS).
def_inverse(diffusion, PS0, PS, "2|s⟩⟨s|−I = H^n (2|0⟩⟨0|−I) H^n, self-inverse") :- !,
    PS0 = ps(N, _, _, _, _, _), N1 is N - 1, numlist(0, N1, Qs),
    foldl([Qa, Sa, Sb]>>apply_ps(h(Qa), Sa, Sb), Qs, PS0, PS1),
    PS1 = ps(_, O, _, _, _, _),
    foldl([Gz, Az0, Az1]>>(lift(Gz, Lz), poly_scale(Lz, -1, NLz), poly_add([[]-1], NLz, Fz), poly_mul(Az0, Fz, Az1)),
          O, [[]-1], Zero),                       % Π (1 - k_b) = [k = 0]
    poly_scale(Zero, -1, NZ), poly_add([[]-1], NZ, NotZero),
    poly_scale(NotZero, 1 rdiv 2, P), add_phase(P, PS1, PS2),
    foldl([Qc, Sc, Sd]>>apply_ps(h(Qc), Sc, Sd), Qs, PS2, PS).

%!  prove(+Circ, +TargetName, -Proved, -Log) is det.
prove(Circ, Name, Ok, Log) :-
    circuit_pathsum(Circ, PS0),
    pathsum_string(PS0, S0),
    format(string(L0), "  circuit path sum: ~w", [S0]),
    (   def_inverse(Name, PS0, PS1, Desc)
    ->  pathsum_string(PS1, S1),
        format(string(L1), "  compose with U† from the definition: ~w", [Desc]),
        format(string(L2), "    ~w", [S1]),
        reduce_pathsum(PS1, PS2, RL),
        maplist([A, B]>>string_concat("  ", A, B), RL, RL2),
        pathsum_string(PS2, S2),
        format(string(L3), "  normal form: ~w", [S2]),
        (   is_identity(PS2)
        ->  Ok = true,
            PS2 = ps(_, _, Ph, _, _, _),
            (   memberchk([]-G, Ph) -> true ; G = 0 ),
            coef_string(G, GS),
            format(string(L4), "  ∎ identity, so the circuit equals U exactly (global phase ~w turn)", [GS])
        ;   Ok = false,
            L4 = "  not proved: the residual above is where circuit and definition differ"
        ),
        append([[L0, L1, L2], RL2, [L3, L4]], Log)
    ;   Ok = false,
        format(string(L1), "  no path-sum definition for ~w", [Name]),
        Log = [L0, L1]
    ).
