-module(poisson).

-export([rpois/1, rpois_big/1, write_csv/1, write_floor_csv/1, normal/2, exp_rand/0]).

-define(M_1_SQRT_2PI, 0.398942280401432677939946059934).
-define(one_7,        0.1428571428571428571).
-define(one_12,       0.0833333333333333333).
-define(one_24,       0.0416666666666666667).
-define(q, [
  0.6931471805599453,
  0.9333736875190459,
  0.9888777961838675,
  0.9984959252914960,
  0.9998292811061389,
  0.9999833164100727,
  0.9999985691438767,
  0.9999998906925558,
  0.9999999924734159,
  0.9999999995283275,
  0.9999999999728814,
  0.9999999999985598,
  0.9999999999999289,
  0.9999999999999968,
  0.9999999999999999,
  1.0000000000000000
]).


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

write_floor_csv(Samples) ->
    {ok, FD} = file:open("/tmp/esample", [write]),
    StringSamples = string:join(
        [io_lib:format("~.9g", [I]) || I <- Samples],
        "\n"
    ),

    io:fwrite(FD, "~s", [StringSamples]),
    file:close(FD).


floor(X) when X < 0 ->
    T = trunc(X),
    case X - T == 0 of
        true -> T;
        false -> T - 1
    end;
floor(X) ->
    trunc(X).


ceiling(X) when X < 0 ->
    trunc(X);
ceiling(X) ->
    T = trunc(X),
    case X - T == 0 of
        true -> T;
        false -> T + 1
    end.

normal(Mean, Sigma) ->
    Rv1 = random:uniform(),
    Rv2 = random:uniform(),
    Rho = math:sqrt(-2 * math:log(1-Rv2)),
    Rho * math:cos(2 * math:pi() * Rv1) * Sigma + Mean.

exp_rand() ->
    U = random:uniform(),
    A = lists:nth(1, ?q),
    [UNew, ANew] = exp_rand_prepare(U, A),
    UMin = random:uniform(),

    exp_rand_sample(UNew, UMin, ANew).

exp_rand_prepare(U, A) when U =< 1.0 ->
    error_logger:info_msg("~p~n", ["prepare"]),
    UNew = U + U,
    ANew = A + lists:nth(1, ?q),
    exp_rand_prepare(UNew, ANew);

exp_rand_prepare(U, A) ->
    [U-1.0, A].

exp_rand_sample(U, UMin, A) ->
    exp_rand_sample(U, UMin, A, 1).


exp_rand_sample(U, UMin, A, Index) ->
    UStar    = random:uniform(),
    CurrentQ = lists:nth(Index, ?q),
    if
        U > CurrentQ ->
            UMinNew = UStar,
            exp_rand_sample(U, UMinNew, A, Index+1);
        true ->
            A + UMin * lists:nth(1, ?q)
    end.












































rpois_big(Mu) ->
    S    = math:sqrt(Mu),
    D    = 6.0 * Mu * Mu,
    BigL = floor(Mu - 1.1484),

    step_n(S, D, Mu, BigL).

step_n(S, D, Mu, BigL) ->
    %error_logger:info_msg("~p~n", ["Step N"]),
    G = Mu + (S * normal(0, 1)),

    if G >= 0.0 ->
        Pois = floor(G),

        if Pois >= BigL -> Pois;
           true         -> step_s(S, D, Mu, Pois)
        end
    end.

step_s(S, D, Mu, Pois) ->
    Fk     = Pois,
    Difmuk = Mu - Fk,
    U      = random:uniform(),

    if (D * U) >= (Difmuk * Difmuk * Difmuk) -> Pois;
       true                                  -> step_p(S, Mu)
    end.

step_p(S, Mu) ->
    error_logger:info_msg("~p~n", ["Step P"]),
    Omega = ?M_1_SQRT_2PI / S,
    % The quantities b1, b2, c3, c2, c1, c0 are for the Hermite
    % approximations to the discrete normal probabilities fk.

    B1 = ?one_24 / Mu,
    B2 = 0.3 * B1 * B1,
    C3 = ?one_7 * B1 * B2,
    C2 = B2 - 15.0 * C3,
    C1 = B1 - 6.0 * B2 + 45.0 * C3,
    C0 = 1.0 - B1 + 3.0 * B2 - 15.0 * C3,
    C  = 0.1069 / Mu,
    200.























