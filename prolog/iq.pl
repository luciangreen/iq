:- encoding(utf8).
:- module(iq, [iq/1, iq/2, iq/3, targets/0]).
/** <module> iq -- quantum algorithms from first principles

Convenience entry points:

  ?- iq(qft).             % 3-qubit QFT, full derivation
  ?- iq(qft, 5).
  ?- iq(toffoli, _, [basis(cx)]).
  ?- targets.
*/
:- reexport(derive).
:- reexport(simplify).
:- reexport(search).
:- reexport(sim).
:- reexport(targets).
:- reexport(output).
:- reexport(poly).
:- reexport(gates).

iq(Name) :- iq(Name, _, []).
iq(Name, N) :- iq(Name, N, []).
iq(Name, N, Opts) :-
    derive(Name, N, Opts, D),
    print_derivation(D, Opts).

targets :-
    forall(target_info(Name, DN, Fixed, Desc),
           (   Fixed == true
           ->  format("  ~w~t~14|~w~n", [Name, Desc])
           ;   format("  ~w N~t~14|~w (default N = ~d)~n", [Name, Desc, DN])
           )).
