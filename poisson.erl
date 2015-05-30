-module(poisson).

-export([
    rpois/1,
    rpois_big/1,
    write_csv/1,
    write_floor_csv/1,
    normal/2,
    exp_rand/0,
    floor/1,
    ceiling/1,
    fsign/2,
    fact/1]).

-define(M_1_SQRT_2PI, 0.398942280401432677939946059934).

-define(a0, -0.5      ).
-define(a1,  0.3333333).
-define(a2, -0.2500068).
-define(a3,  0.2000118).
-define(a4, -0.1661269).
-define(a5,  0.1421878).
-define(a6, -0.1384794).
-define(a7,  0.1250060).

-define(one_7,        0.1428571428571428571).
-define(one_12,       0.0833333333333333333).
-define(one_24,       0.0416666666666666667).
-define(Q0,           0.6931471805599453).
-define(Q, [
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
-define(fact, [1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0]).


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

fsign(X, Y) when Y >= 0 ->
    erlang:abs(X);
fsign(X, _) ->
    -erlang:abs(X).


fact(N) -> fact(N,1).

fact(0,Acc) -> Acc;
fact(N,Acc) when N > 0 -> fact(N-1,N*Acc).


normal(Mean, Sigma) ->
    Rv1 = random:uniform(),
    Rv2 = random:uniform(),
    Rho = math:sqrt(-2 * math:log(1-Rv2)),
    Rho * math:cos(2 * math:pi() * Rv1) * Sigma + Mean.

exp_rand() ->
    A = 0.0,
    U = random:uniform(),

    [UNew, ANew] = exp_rand_prepare(U, A),

    if
        UNew < ?Q0 ->
            UNew + ANew;
        true ->
            UMin = random:uniform(),
            exp_rand_sample(UNew, UMin, ANew, ?Q)
    end.


exp_rand_prepare(U, A) ->
    UNew = U + U,

    if UNew > 1.0 -> [UNew - 1.0, A];
        true      -> exp_rand_prepare(UNew, A + ?Q0)
    end.


exp_rand_sample(U, UMin, A, [Q|QRest]) when U > Q ->
    UStar = random:uniform(),

    if
        UMin > UStar -> exp_rand_sample(U, UStar, A, QRest);
        true         -> exp_rand_sample(U, UMin,  A, QRest)
    end;


exp_rand_sample(_, UMin, A, _) ->
    A + UMin * ?Q0.



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
           true         -> step_s(S, D, G, Mu, Pois)
        end
    end.

step_s(S, D, G, Mu, Pois) ->
    Fk     = Pois,
    Difmuk = Mu - Fk,
    U      = random:uniform(),

    if (D * U) >= (Difmuk * Difmuk * Difmuk) -> Pois;
       true                                  -> step_p(Pois, S, G, Mu, Difmuk, Fk, U)
    end.

step_p(Pois, S, G, Mu, Difmuk, Fk, U) ->
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

    if
        G >= 0.0 -> step_f1(Pois, Mu, Fk, Difmuk, nan, U, S, Omega, C, C0, C1, C2, C3);
        true     -> step_e(Pois, Mu, Fk, Difmuk, S, Omega, C, C0, C1, C2, C3)
    end.

step_e(Pois, Mu, Fk, Difmuk, S, Omega, C, C0, C1, C2, C3) ->
    E = exp_rand(),
    U = 2 * random:uniform() - 1,
    T = fsign(E, U),
    step_f0(Pois, Mu, Fk, Difmuk, E, U, S, T, Omega, C, C0, C1, C2, C3).


step_f0(_Pois, Mu, _Fk, _Difmuk, E, U, S, T, Omega, C, C0, C1, C2, C3) when -0.6744 < T ->
    Pois   = floor(Mu + (S * T)),
    Fk     = Pois,
    Difmuk = Mu - Fk,
    step_f1(Pois, Mu, Fk, Difmuk, E, U, S, Omega, C, C0, C1, C2, C3);


step_f0(Pois, _, _, _, _, _, _, _, _, _, _, _, _, _) ->
    Pois.


step_f1(Pois, Mu, Fk, Difmuk, E, U, S, Omega, C, C0, C1, C2, C3) ->
    {Px, Py} = if
        Pois < 10 ->
            {-Mu, ( math:pow(Mu, Pois) / lists:nth(Pois, ?fact) )};
        true ->
            Del0  = ?one_12 / Fk,
            Del1  = Del0 * (1.0 - 4.8 * Del0 * Del0),
            V     = Difmuk / Fk,
            FabsV = erlang:abs(V),

            if
                 FabsV =< 0.25 ->
                    error_logger:info_msg("STEP SUB", []),
                    PxTemp = Fk * V * V * (
                        (
                            (
                                (
                                    (
                                        (
                                            (
                                                ?a7 * V + ?a6
                                            ) * V + ?a5
                                        ) * V + ?a4
                                    ) * V + ?a3
                                ) * V + ?a2
                            ) * V + ?a1
                        ) * V + ?a0
                    ) - Del1,
                    {PxTemp, 0.0};
                true ->
                    {
                        Fk * math:log(1.0 + V) - Difmuk - Del1,
                        ?M_1_SQRT_2PI / math:sqrt(Fk)
                    }
            end
    end,


    X0 = (0.5 - Difmuk) / S,
    X1 = X0 * X0,
    Fx = -0.5 * X1,
    Fy = Omega * (((C3 * X1 + C2) * X1 + C1) * X1 + C0),

    case E of
        nan ->
            error_logger:info_msg("~p~n", [{Fx, Fy, U, Px, Py}]),
            case Fy - U * Fy =< Py * math:exp(Px - Fx) of
                true -> Pois;
                _    -> step_e(Pois, Mu, Fk, Difmuk, S, Omega, C, C0, C1, C2, C3)
            end;
        _ ->
            case C * erlang:abs(U) =< Py * math:exp(Px + E) - Fy * math:exp(Fx + E) of
                true -> Pois;
                _    -> step_e(Pois, Mu, Fk, Difmuk, S, Omega, C, C0, C1, C2, C3)
            end
    end.






