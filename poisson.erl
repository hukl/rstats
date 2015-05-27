-module(poisson).

-export([rpois/1, write_csv/1]).


rpois(Mu) ->
    M  = erlang:max(1, trunc(Mu)),
    Q  = P0 = P = math:exp(-Mu),
    U  = random:uniform(),
    PP = [],
    if U =< P0 -> 0;
       true    -> random_poisson(1, U, M, Mu, P, PP, Q)
    end.


lookup_poisson(M, PP) ->
    U  = random:uniform(),
    MinStart = case U =< 0.458 of
        true -> 1;
        _    -> erlang:min(35, M)
    end,

    FilterFun = fun(Element) ->
        U =< lists:nth(Element, PP)
    end,

    Result = lists:filter(FilterFun, lists:seq(MinStart, 35)),

    case Result of
        [Head|_] -> Head;
        []       -> lookup_poisson(M, PP)
    end.

random_poisson(36, _, M, _, _, PP, _) ->
    lookup_poisson(M, PP);

random_poisson(K, U, M, Mu, P, PP, Q) ->
    PNew  = P * Mu / K,
    QNew  = Q + PNew,

    PPNew = lists:append(PP, [QNew]),
    if U =< QNew -> K;
       true       -> random_poisson(K+1, U, M, Mu, PNew, PPNew, QNew)
    end.

write_csv(Samples) ->
    {ok, FD} = file:open("/tmp/esample", [write]),
    StringSamples = string:join(
        [integer_to_list(I) || I <- Samples],
        "\n"
    ),

    io:fwrite(FD, "~s", [StringSamples]),
    file:close(FD).
