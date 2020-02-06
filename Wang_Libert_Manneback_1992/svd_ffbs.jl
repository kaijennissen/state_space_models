#------------------------------------------------------------------------------
# Julia Code for a SVD-Kalman-Filter
# Wang, L., Libert, G. & Manneback, P. (1992).
# Kalman Filter Algorithm Based on Singualar Value Decomposition.
#------------------------------------------------------------------------------
using LinearAlgebra
using CSV
using Plots
using Statistics
using Distributions



function kalman_filter(Y, G, F, W, V, C0, m0)
    T = size(Y, 1);
    q1 = size(V, 1);
    q2 = size(W, 1);

    store_theta = zeros(q2, T);

    T = size(y, 1);
    m = m0;
    C = C0;

    for t in 1:T

        # pred state
        a = G*m
        R = G*C*G'+W

        # pred observation
        f = F*a
        Q = F*R*F'+V

        # filtered
        e = Y[t:t,]-f
        m = a + R*F'*inv(Q)*e
        C = R-R*F'*inv(Q)*F*R

        store_theta[:,t] = m
    end
    return store_theta
end

function svd_kalman(Y, G, F, W, V, C0, m0)
    T = size(Y, 1);
    q1 = size(V, 1);
    q2 = size(W, 1);

    store_m = zeros(q2, T);

    m = m0;
    C = C0;
    xm = m0;
    xC = C0;

    USVT = svd(C);
    U_plus = Matrix(USVT.U);
    D_plus = Diagonal(USVT.S)^0.5;

    for t in 1:T

        # predict
        a = G*m

        # update U and D
        W_sq = W^0.5
        USVT = svd(Matrix([D_plus*U_plus'*G' ; W_sq]))
        U = USVT.V
        D = Diagonal(USVT.S)

        L = Matrix(cholesky(inv(V)).L)

        USVT = svd(Matrix([L'*F*U; D^-1]))
        V_star = USVT.V
        D_star = Diagonal(USVT.S)

        U_plus = U*V_star
        D_plus = D_star^-1

        # compute gain K
        K = U_plus*(D_plus^2)*U_plus'*F'*L*L'

        # update estimate x
        m = a + K*(Y[t:t]-F*a)

        store_m[:,t] = m;

    end
    return store_m
end


function svd_ffbs(Y, G, F, W, V, C0, m0)
    epss = eps(Float64)^0.4
    T = size(Y, 1);
    q1 = size(V, 1);
    q2 = size(W, 1);

    store_xa = zeros(q2, T);
    store_xm = zeros(q2, T+1);
    store_a = zeros(q2, T);
    store_m = zeros(q2, T+1);
    store_D_plus = zeros(q2, q2, T+1);
    store_U_plus = zeros(q2, q2, T+1);
    store_D = zeros(q2, q2, T);
    store_U = zeros(q2, q2, T);
    store_theta = zeros(q2, T+1)

    m = m0;
    C = C0;
    xm = m0;
    xC = C0;
    # predict
    a = G*m

    tmp = svd(C);
    U_plus = Matrix(tmp.Vt);
    D_plus = Diagonal(tmp.S .^0.5);
    #D_inv = Diagonal(tmp.S .^-0.5);

    store_xm[:,1] = m;
    store_a[:,1] = a;
    store_m[:,1] = m;
    store_D_plus[:,:,1] = D_plus;
    store_U_plus[:,:,1] = U_plus;

    for t in 1:T

        tmp = svd(V);
        Uv = tmp.Vt';
        Dv = tmp.S.^0.5
        Dv_inv = 1 ./ Dv;
        Dv_inv[findall(x -> x == Inf, Dv_inv)] .= 0
        sqrtVinv = Dv_inv * tmp.Vt

        tmp = svd(W)
        sqrtW = Diagonal(tmp.S)^0.5 * tmp.Vt

        # prior
        a = G*m
        tmp = svd(Matrix([D_plus*U_plus'*G'; sqrtW]))
        U = tmp.Vt'
        D = Diagonal(tmp.S)
        D_inv = Diagonal(tmp.S .^-1)
        D_inv[findall(x -> x == Inf, D_inv)] .= 0

        # one-step forecast
        f = F * a

        # posterior
        tmp = svd(Matrix([sqrtVinv*F*U; D_inv]));
        U_plus = U*tmp.Vt';   #U*V_star;
        D_plus = Diagonal(tmp.S .^-1);
        D_plus[findall(x -> x == Inf, D_plus)] .= 0;
        #V_star = tmp.Vt;
        #D_star = Diagonal(tmp.S);
        #D_star_inv = Diagonal(tmp.S .^-1);
        #D_star_inv[findall(x -> x == Inf, D_star_inv)] .= 0

        # compute gain K
        K = U_plus*(D_plus*(D_plus*(U_plus'*(F'*(sqrtVinv'*sqrtVinv)))));

        # update estimate
        m = a + K*(Y[t:t]-F*a);

        # pred state
        xa = G*xm
        xR = G*xC*G' + W
        xK = xR*F'*inv(F*xR*F'+V)
        xm = xa + xR*F'*inv(F*xR*F'+V)*(Y[t:t,]-F*xa)
        xC = xR-xR*F'*inv(F*xR*F'+V)*F*xR


        store_xa[:,t] = xa;
        store_xm[:,t+1] = xm;

        store_a[:,t] = a;
        store_m[:,t+1] = m;
        store_D_plus[:,:,t+1] = D_plus;
        store_U_plus[:,:,t+1] = U_plus;
        store_D[:,:,t] = D;
        store_U[:,:,t] = U;
    end

    # m = store_m[:,T+1];
    # D_plus = store_D_plus[:,:,T+1];
    # U_plus = store_U_plus[:,:,T+1];
    #
    # h = m;
    #
    # theta = h + U_plus*(D_plus*randn(q2))
    # store_theta[:,T+1] = theta;
    #
    # for t in collect(T:-1:1)
    #
    #     a = store_a[:,t];
    #     m = store_m[:,t];
    #     D_plus = store_D_plus[:,:,t];
    #     U_plus = store_U_plus[:,:,t];
    #     D = store_D[:,:,t];
    #     U = store_U[:,:,t];
    #
    #     tmp = svd(W)
    #     DW = tmp.S .^ 0.5
    #     DW = maximum([DW epss*ones(13,1)], dims=2)
    #     DW_inv = Diagonal(1 ./ DW)
    #     sqrtWinv = tmp.Vt'*DW_inv^2*tmp.Vt
    #
    #     L_star = DW_inv*tmp.Vt
    #     tmp = svd(Matrix([L_star'*G*U_plus; D_plus^-1]))
    #     V_tri = tmp.Vt
    #
    #     U_sq = U_plus*V_tri'
    #     D_sq = tmp.S .^ -1
    #
    #     D_sq[findall(a-> a == Inf, D_sq)] .= 0
    #     D_sq = Diagonal(D_sq)
    #
    #     h = m + U_sq*D_sq^2*U_sq'*G'sqrtWinv*(theta-a)
    #     theta = h + U_sq*D_sq*randn(q2)
    #
    #     store_theta[:,t] = theta
    #
    # end

    return store_m', store_xm'
    #return store_theta
end

function local_trend_seasonal(y)
    Y = broadcast(log, y)
    q = 13
    F = [1 0 1 zeros(1, 10)];
    G = zeros(q, q);
    G[1, 1] = 1;
    G[1, 2] = 1;
    G[2, 2] = 1;
    G[3, 3:end] = -1*ones(11);
    G[4:end, 3:end-1] = 1.0I(10);

    V = 1.0/22814.31*ones(1, 1);
    W = Diagonal([1.0/2472.64; 1.0/55947.63; 1.0/4309.207, zeros(10,1)]); #zeros(q, q);
    #W[1, 1] = 1.0/2472.64;
    #W[2, 2] = 1.0/55947.63;
    #W[3, 3] = 1.0/4309.207;

    m0 = zeros(q,1);
    C0 = 1e7I(q);
    return Y, G, F, W, V, C0, m0
end


# Example 1
data_raw = CSV.read("./bayesian_inference/AirPassengers.csv", header = 0);
y = map(x->parse(Float64,x), data_raw[2:end, 2]);
plot(y);

Y, G, F, W, V, C0, m0 = local_trend_seasonal(y);

# Conventional Kalman Filter
store_theta = kalman_filter(Y, G, F, W, V, C0, m0)
plot(Y, legend=false)
plot!(store_theta[1,:])

plot(Y)
plot!((F*store_theta)')


# Singular Value Algorithm
store_theta = svd_kalman(Y, G, F, W, V, C0, m0)

plot(Y, legend=false)
plot!(store_theta[1,:])

plot(Y)
plot!((F*store_theta)')


# Singular Value Algorithm - FFBS
theta_svd, theta = svd_ffbs(Y, G, F, W, V, C0, m0)

plot(Y, legend=false)
plot!(theta[1,:])

plot(Y)
plot!((F*theta)')
