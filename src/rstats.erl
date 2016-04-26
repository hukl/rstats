-module(rstats).

-export([
    rpois/1,
    write_csv/2,
    write_float_csv/2,
    normal/2,
    exp_rand/0,
    rexp/1,
    walker_lookup_table/1,
    walker_choice/1,
    ewa/2,
    ewa_next_state/3,
    floor/1,
    ceiling/1,
    fsign/2,
    fact/1,
    dnorm/1,
    dnorm/2,
    dnorm/3,
    rtruncnorm/5]).

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

-define(M_LN2, 0.693147180559945309417232121458).

-define(DBL_MANT_DIG, 16).

-define(DBL_MIN_EXP, -999).


-define(t1, 0.15).
-define(t2, 2.18).
-define(t3, 0.725).
-define(t4, 0.45).

% REFERENCE
%
% Ahrens, J.H. and Dieter, U. (1982).
% Computer generation of Poisson deviates
% from modified normal distributions.
% ACM Trans. Math. Software 8, 163-179.
%
% Re-Implemented from R's nmath/rpois.c C Code

rpois(Mu) when Mu < 10.0 ->
    M  = erlang:max(1, trunc(Mu)),
    Q  = P0 = P = math:exp(-Mu),
    U  = random:uniform(),
    PP = [],
    if U =< P0 -> 0;
       true    -> step_c(1, U, M, Mu, P, PP, Q)
    end;

rpois(Mu) ->
    S    = math:sqrt(Mu),
    D    = 6.0 * Mu * Mu,
    BigL = floor(Mu - 1.1484),

    step_n(S, D, Mu, BigL).


step_u(M, PP) ->
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
        []       -> step_u(M, PP)
    end.


step_c(36, _, M, _, _, PP, _) ->
    step_u(M, PP);

step_c(K, U, M, Mu, P, PP, Q) ->
    PNew  = P * Mu / K,
    QNew  = Q + PNew,

    PPNew = lists:append(PP, [QNew]),
    if U =< QNew -> K;
       true       -> step_c(K+1, U, M, Mu, PNew, PPNew, QNew)
    end.


step_n(S, D, Mu, BigL) ->
    G = Mu + (S * normal(0, 1)),

    Pois = if
        G >= 0.0 -> floor(G);
        true     -> -1
    end,

    if
        Pois >= BigL -> Pois;
        true         -> step_s(S, D, G, Mu, Pois)
    end.


step_s(S, D, G, Mu, Pois) ->
    Fk     = Pois,
    Difmuk = Mu - Fk,
    U      = random:uniform(),

    if (D * U) >= (Difmuk * Difmuk * Difmuk) -> Pois;
       true                                  -> step_p(Pois, S, G, Mu, Difmuk, Fk, U)
    end.


step_p(Pois, S, G, Mu, Difmuk, Fk, U) ->
    Omega = ?M_1_SQRT_2PI / S,
    B1    = ?one_24 / Mu,
    B2    = 0.3 * B1 * B1,
    C3    = ?one_7 * B1 * B2,
    C2    = B2 - 15.0 * C3,
    C1    = B1 - 6.0 * B2 + 45.0 * C3,
    C0    = 1.0 - B1 + 3.0 * B2 - 15.0 * C3,
    C     = 0.1069 / Mu,

    if
        G >= 0.0 -> step_f1(Pois, Mu, Fk, Difmuk, nan, U, S, Omega, C, C0, C1, C2, C3);
        true     -> step_e(Pois, Mu, Fk, Difmuk, S, Omega, C, C0, C1, C2, C3)
    end.


step_e(Pois, Mu, Fk, Difmuk, S, Omega, C, C0, C1, C2, C3) ->
    E = exp_rand(),
    U = (2 * random:uniform()) - 1.0,
    T = 1.8 + fsign(E, U),
    step_f0(Pois, Mu, Fk, Difmuk, E, U, S, T, Omega, C, C0, C1, C2, C3).


step_f0(_Pois, Mu, _Fk, _Difmuk, E, U, S, T, Omega, C, C0, C1, C2, C3) when -0.6744 < T ->
    Pois   = floor(Mu + S * T),
    Fk     = Pois,
    Difmuk = Mu - Fk,
    step_f1(Pois, Mu, Fk, Difmuk, E, U, S, Omega, C, C0, C1, C2, C3);

step_f0(Pois, Mu, Fk, Difmuk, _E, _U, S, _T, Omega, C, C0, C1, C2, C3) ->
    step_e(Pois, Mu, Fk, Difmuk, S, Omega, C, C0, C1, C2, C3).


step_f1(Pois, Mu, Fk, Difmuk, E, U, S, Omega, C, C0, C1, C2, C3) ->
    {Px, Py} = if
        Pois < 10 ->
            {-Mu, ( math:pow(Mu, Pois) / lists:nth((Pois+1), ?fact) )};
        true ->
            Del0   = ?one_12 / Fk,
            Del1   = Del0 * (1.0 - 4.8 * Del0 * Del0),
            V      = Difmuk / Fk,
            FabsV  = erlang:abs(V),
            PyTemp = ?M_1_SQRT_2PI / math:sqrt(Fk),

            if
                 FabsV =< 0.25 ->
                    PxTemp = Fk * V * V * (
                        ((((((?a7 * V + ?a6) * V + ?a5) * V + ?a4) * V + ?a3) * V + ?a2) * V + ?a1) * V + ?a0
                    ) - Del1,
                    {PxTemp, PyTemp};
                true ->
                    PxTemp = Fk * math:log(1.0 + V) - Difmuk - Del1,
                    {PxTemp, PyTemp}
            end
    end,


    X0 = (0.5 - Difmuk) / S,
    X1 = X0 * X0,
    Fx = -0.5 * X1,
    Fy = Omega * (((C3 * X1 + C2) * X1 + C1) * X1 + C0),

    case E of
        nan ->
            case (Fy - U * Fy) =< (Py * math:exp(Px - Fx)) of
                true -> Pois;
                _    -> step_e(Pois, Mu, Fk, Difmuk, S, Omega, C, C0, C1, C2, C3)
            end;
        _ ->
            case (C * erlang:abs(U)) =< (Py * math:exp(Px + E) - Fy * math:exp(Fx + E)) of
                true -> Pois;
                _    -> step_e(Pois, Mu, Fk, Difmuk, S, Omega, C, C0, C1, C2, C3)
            end
    end.


% REFERENCE
%
% Ahrens, J.H. and Dieter, U. (1972).
% Computer methods for sampling from the exponential and
% normal distributions.
% Comm. ACM, 15, 873-882.
%
% Re-Implemented from R's nmath/sexp.c C Code

% Random variates from the standard exponential distribution
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


% Random variates from the exponential distribution
rexp(Scale) when 0 < Scale ->
    Scale * exp_rand().


% Normal Distribution / Box-Muller
normal(Mean, Sigma) ->
    Rv1 = random:uniform(),
    Rv2 = random:uniform(),
    Rho = math:sqrt(-2 * math:log(1-Rv2)),
    Rho * math:cos(2 * math:pi() * Rv1) * Sigma + Mean.



% Compute the density of the normal distribution.
dnorm(X) ->
    % (x, mean, sigma)
    dnorm(X, 0, 1).

dnorm(X, Mean) ->
    dnorm(X, Mean, 1).

dnorm(_X, _Mean, Sigma) when Sigma =< 0 ->
    {error, not_a_number};

dnorm(X, Mean, Sigma) ->
    X1 = (X - Mean) / Sigma,

    X2 = erlang:abs(X1),

    do_dnorm(X2, Sigma).

do_dnorm(X, Sigma) when X < 5 ->
    ?M_1_SQRT_2PI * math:exp(-0.5 * X * X) / Sigma;

do_dnorm(X, Sigma) ->
    case (X > math:sqrt(-2 * ?M_LN2 * (?DBL_MIN_EXP + 1 - ?DBL_MANT_DIG))) of
        true ->
            0.0;
        false ->
            X1 = ldexp(
                trunc(ldexp(X, 16)), -16
            ),

            X2 = X - X1,
            ?M_1_SQRT_2PI / Sigma * ( math:exp(-0.5 * X1 * X1) * math:exp( (-0.5 * X2 - X1) * X2 ) )
    end.





ers_a_inf(A) ->
    Ainv = 1.0 / A,
    ers_a_inf(Ainv, A).

ers_a_inf(Ainv, A) ->
    X   = rexp(Ainv) + A,
    Rho = math:exp(-0.5 * math:pow((X - A), 2)),
    R   = crypto:rand_uniform(0, 1),

    case R > Rho of
        true -> ers_a_inf(Ainv, A);
        _    -> X
    end.


ers_a_b(A, B) ->
    Ainv = 1.0 / A,
    ers_a_b(Ainv, A, B).

ers_a_b(Ainv, A, B) ->
    X   = rexp(Ainv) + A,
    Rho = math:exp(-0.5 * math:pow((X - A), 2)),
    R   = crypto:rand_uniform(0, 1),

    case (R > Rho) orelse (X > B) of
        true -> ers_a_b(Ainv, A, B);
        _    -> X
    end.



nrs_a_inf(A) ->
    X = normal(0, 1),

    case X < A of
        true -> nrs_a_inf(A);
        _    -> X
    end.


nrs_a_b(A, B) ->
    X = normal(0, 1),

    case (X < A) orelse (X > B) of
        true -> nrs_a_b(A, B);
        _    -> X
    end.


hnrs_a_b(A, B) ->
    X1 = normal(0, 1),
    X2 = erlang:abs(X1),

    case (X2 < A) orelse (X2 > B) of
        true -> hnrs_a_b(A, B);
        _    -> X2
    end.


urs_a_b(A, B) ->
    Phi_a = dnorm(A, 0.0, 1.0),

    Ub = case (A < 0) andalso (B > 0) of
        true -> ?M_1_SQRT_2PI;
        _    -> Phi_a
    end,

    do_urs_a_b(A, B, Ub, Phi_a).

do_urs_a_b(A, B, Ub, Phi_a) ->
    X = crypto:rand_uniform(A, B),

    case crypto:rand_uniform(0, 1) * Ub > dnorm(X, 0, 1) of
        true -> do_urs_a_b(A, B, Ub, Phi_a);
        _    -> X
    end.


r_lefttruncnorm(A, Mean, Sd) ->
    Alpha = (A - Mean) / Sd,

    case Alpha < ?t4 of
        true -> Mean + Sd * nrs_a_inf(Alpha);
        _    -> Mean + Sd * ers_a_inf(Alpha)
    end.


r_righttruncnorm(B, Mean, Sd) ->
    Beta = (B - Mean) / Sd,
    Mean - Sd * r_lefttruncnorm(-Beta, 0.0, 1.0).



% NOTE: if something is not right check if runif is really crypto:rand_uniform
r_truncnorm(A, B, Mean, Sd) ->
    Alpha = (A - Mean) / Sd,
    Beta = (B - Mean) / Sd,

    Phi_a = dnorm(Alpha, 0.0, 1.0),
    Phi_b = dnorm(Beta, 0.0, 1.0),

    if
        Beta =< Alpha ->
            {error, na_real};
        (Alpha =< 0) andalso (0 =< Beta) ->
            if
                (Phi_a =< ?t1) orelse (Phi_b =< ?t1) ->
                    Mean + Sd * nrs_a_b(Alpha, Beta);
                true ->
                    Mean + Sd * urs_a_b(Alpha, Beta)
            end;
        Alpha > 0 ->
            if
                (Phi_a / Phi_b =< ?t2) ->
                    Mean + Sd * urs_a_b(Alpha, Beta);
                true ->
                    if
                        Alpha < ?t3 ->
                            Mean + Sd * hnrs_a_b(Alpha, Beta);
                        true ->
                            Mean + Sd * ers_a_b(Alpha, Beta)
                    end
            end;
        true ->
            if
                Phi_b / Phi_a =< ?t2 ->
                    Mean - Sd * urs_a_b(-Beta, -Alpha);
                true ->
                    if
                        Beta > -(?t3) ->
                            Mean - Sd * hnrs_a_b(-Beta, -Alpha);
                        true ->
                            Mean - Sd * ers_a_b(-Beta, -Alpha)
                    end
            end
    end.


% Truncated normal distribution
rtruncnorm(N, infinity, infinity, Mean, Sd) ->
    [normal(Mean, Sd) || _ <- lists:seq(1, N)];

rtruncnorm(N, infinity, B, Mean, Sd) ->
    [r_righttruncnorm(B, Mean, Sd) || _ <- lists:seq(1, N)];

rtruncnorm(N, A, infinity, Mean, Sd) ->
    [r_lefttruncnorm(A, Mean, Sd)  || _ <- lists:seq(1, N)];

rtruncnorm(N, A, B, Mean, Sd) ->
    [r_truncnorm(A, B, Mean, Sd)  || _ <- lists:seq(1, N)].


% Walker alias method - efficient random selection with defined probabilities.
% <http://en.wikipedia.org/wiki/Alias_method>
%
% > ProbabilityValueList = [{10, a}, {20, b}, {30, c}, {40, d}],
% > WalkerVectors = walker_alias:build(ProbabilityValueList),
% > RandomValue = walker_alias:choice(WalkerVectors).
%
% Probability values may be any integers
% RandomValue will contain 'd' in 40% cases, 'c' in 30% and so on.
% WalkerVectors is reusable struct - you can call choice/1 on it for many
% times.
%
% Ported from Python implementation:
% <http://code.activestate.com/recipes/576564-walkers-alias-method-for-random-objects-with-diffe/>
% Other implementations:
% Common Lisp: http://prxq.wordpress.com/2006/04/23/more-on-the-alias-method/
% Ruby: https://github.com/cantino/walker_method
% Created : 16 May 2013 by Sergey Prokhorov <me@seriyps.ru>
%
walker_lookup_table(WeightedList) ->
    {Keys, Weights} = lists:unzip(WeightedList),
    walker_lookup_table(Keys, Weights).


walker_lookup_table(Keys, Weights) ->
    N             = length(Weights),
    Sumw          = lists:sum(Weights),
    Prob          = [W * N / Sumw || W <- Weights],
    {Short, Long} = split_short_long(Prob),
    {Inx, Prob2}  = build_vectors(Short, Long, array:new(N), array:from_list(Prob)),

    {walker_vectors, N, array:from_list(Keys), Inx, Prob2}.


split_short_long(Probilities) ->
    SplitFun = fun(Probability, {Sh, Lo, Index}) when Probability < 1 ->
                       {[Index | Sh], Lo, Index + 1};
                  (Probability, {Sh, Lo, Index}) when Probability >= 1 ->
                       {Sh, [Index | Lo], Index + 1}
               end,

    {Short, Long, _} = lists:foldl(SplitFun, {[], [], 0}, Probilities),
    {Short, Long}.


build_vectors(S, L, Inx, Prob) when (S == []) orelse (L == []) ->
    {Inx, Prob};

build_vectors([J | Short], [K | Long], Inx, Prob) ->
    NewInx   = array:set(J, K, Inx),
    ProbJ    = array:get(J, Prob),
    ProbK    = array:get(K, Prob),
    NewProbK = ProbK - (1 - ProbJ),
    NewProb  = array:set(K, NewProbK, Prob),

    {NewShort, NewLong} = case NewProbK < 1 of
        true  -> {[K | Short], Long};
        false -> {Short, [K | Long]}
    end,

    build_vectors(NewShort, NewLong, NewInx, NewProb).


walker_choice({walker_vectors, N, Keys, Inx, Prob}) ->
    J = crypto:rand_uniform(0, N),
    U = crypto:rand_uniform(0, 100) / 100,
    Idx = case U =< array:get(J, Prob) of
              true -> J;
              false -> array:get(J, Inx)
          end,
    array:get(Idx, Keys).



%% Exponentially Weighted Moving Average
%% https://en.wikipedia.org/wiki/Moving_average#Exponential_moving_average
ewa(_, []) -> 0;
ewa(Alpha, [Head|Tail]) -> ewa(Alpha, (1 - Alpha), Tail, Head).

ewa(_, _, [], Acc) -> Acc;
ewa(Alpha, AlphaI, [Head|Tail], Acc) ->
    NewAcc = Alpha * Head + AlphaI * Acc,
    ewa(Alpha, AlphaI, Tail, NewAcc).

ewa_next_state(Value, Alpha, State) ->
    Alpha * Value + (1 - Alpha) * State.

% Helpers missing in Erlangs Standard Library

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


ldexp(X, Exponent) ->
    X * math:pow(2, Exponent).


% Helpers for verifying / comparing results in R

write_csv(Path, Samples) ->
    {ok, FD} = file:open(Path, [write]),
    StringSamples = string:join(
        [integer_to_list(I) || I <- Samples],
        "\n"
    ),

    io:fwrite(FD, "~s", [StringSamples]),
    file:close(FD).

write_float_csv(Path, Samples) ->
    {ok, FD} = file:open(Path, [write]),
    StringSamples = string:join(
        [io_lib:format("~.9g", [I]) || I <- Samples],
        "\n"
    ),

    io:fwrite(FD, "~s", [StringSamples]),
    file:close(FD).

