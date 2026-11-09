:- encoding(utf8).
/*  Command line:

    swipl prolog/cli.pl qft 4
    swipl prolog/cli.pl qft 8 --approx 0.05 --quiet
    swipl prolog/cli.pl toffoli --basis cx --qasm
    swipl prolog/cli.pl iqft 3 --no-swaps
    swipl prolog/cli.pl targets
*/
:- use_module(library(main)).
:- use_module(iq).
:- initialization(main, main).

main(Argv) :-
    (   Argv = [] ; Argv = [help|_] ; Argv = ['--help'|_] )
    ->  usage
    ;   Argv = [targets|_]
    ->  targets
    ;   Argv = [NameA|Rest],
        atom_string(NameA, NS), atom_string(Name, NS),
        parse(Rest, N, Opts),
        catch(iq(Name, N, Opts), E, report(E))
    ).

parse([], _, []).
parse(['--basis', B|T], N, [basis(B)|O]) :- !, parse(T, N, O).
parse(['--approx', E|T], N, [approx(V)|O]) :- !, atom_number(E, V), parse(T, N, O).
parse(['--no-swaps'|T], N, [swaps(false)|O]) :- !, parse(T, N, O).
parse(['--no-simplify'|T], N, [simplify(false)|O]) :- !, parse(T, N, O).
parse(['--qasm'|T], N, [qasm(true)|O]) :- !, parse(T, N, O).
parse(['--quiet'|T], N, [quiet(true)|O]) :- !, parse(T, N, O).
parse(['--no-verify'|T], N, [verify(false)|O]) :- !, parse(T, N, O).
parse(['--verify'|T], N, [verify(true)|O]) :- !, parse(T, N, O).
parse(['--depth', D|T], N, [depth(V)|O]) :- !, atom_number(D, V), parse(T, N, O).
parse([A|T], N, O) :- atom_number(A, N), parse(T, _, O).

usage :-
    format("iq -- derive quantum circuits from first principles~n~n"),
    format("usage: swipl prolog/cli.pl TARGET [N] [options]~n"),
    format("  --basis cx      lower controlled phases to CNOT + phase~n"),
    format("  --approx EPS    drop rotations smaller than EPS radians~n"),
    format("  --no-swaps      leave output qubits permuted instead of SWAPs~n"),
    format("  --no-simplify   show the raw derived circuit~n"),
    format("  --qasm          also print OpenQASM 3~n"),
    format("  --quiet         circuit only, no derivation trace~n"),
    format("  --depth D       search depth for unstructured targets (default 4)~n~n"),
    format("targets:~n"),
    targets.

report(synthesis_error(M)) :- !, format(user_error, "iq: ~w~n", [M]), halt(1).
report(E) :- print_message(error, E), halt(1).
