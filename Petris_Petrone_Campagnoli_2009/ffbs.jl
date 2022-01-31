#------------------------------------------------------------------------------
# Julia Code for a Singular Value Based
# Forward-Filtering Backward-Sampling Algorithm based on:
#
# Petris, G., Petrone, S. & Campagnoli, P. (2009).
# Dynamic Linear Models with R.
#
# Wang, L.,  Libert, G. & Manneback, P. (1992).
# Kalman Filter Algorithm based on Singular Value Decomposition.
# Proceedings of the 31st IEEE Conference on Decision and Control
#
# Kulikova, M. & Tsyganova, J. (2017).
# Improved Discrete-Time Kalman Filtering within Singular Value Decomposition
# IET Control Theory & Applications, 11 (15):2412
#------------------------------------------------------------------------------

# State Space Equation
# y_t = X_t * beta_t + F_t * x_t + vega_t  vega_t ~ N(0, V_t)
# eta_t = Z_t * gamma_t + G_t * eta_t-1 + omega_t  omega_t ~ N(0, W_t)

using LinearAlgebra
using CSV
using Plots
using Statistics
using Distributions
using Random


include("bayesian_utils.jl")
add_dim(x::Array) = reshape(x, (size(x)...,1))

# SVD Forward-Backward-Algorithm
function svd_forward_backward(Y, G, F, W, V, C0, m0)
    epss = eps(Float64)^0.4
    T = size(Y, 1);
    q1 = size(V, 1);
    q2 = size(W, 1);

    store_a = zeros(q2, T);
    store_m = zeros(q2, T+1);
    store_D_plus = zeros(q2, q2, T+1);
    store_U_plus = zeros(q2, q2, T+1);
    store_theta = zeros(q2, T+1)

    m = m0;
    C = C0;

    tmp = svd(C);
    U_plus = Matrix(tmp.Vt);
    D_plus = Diagonal(tmp.S .^0.5);

    tmp = svd(V);
    Uv = tmp.Vt';
    Dv = tmp.S.^0.5
    Dv_inv = 1 ./ Dv;
    Dv_inv[findall(x -> x == Inf, Dv_inv)] .= 0
    sqrtVinv = Dv_inv * tmp.Vt

    tmp = svd(W);
    Uw = tmp.Vt';
    Dw = tmp.S .^0.5
    sqrtW = Diagonal(Dw)
    Dw_inv = 1 ./ Dw ;
    Dw_inv[findall(x -> x == Inf, Dw_inv)] .= 0
    sqrtWinv = Diagonal(Dw_inv) * tmp.Vt

    store_m[:,1] = m;
    store_D_plus[:,:,1] = D_plus;
    store_U_plus[:,:,1] = U_plus;
	
    # forward 
    for t in 1:T

        # prior
        a = G*m
        tmp = svd(Matrix([D_plus*U_plus'*G'; sqrtW]))
        U = tmp.Vt'
        D = Diagonal(tmp.S)
        D_inv = Diagonal(tmp.S .^-1)
        D_inv[findall(x -> x == Inf, D_inv)] .= 0

        # posterior
        tmp = svd(Matrix([sqrtVinv*F*U; D_inv]));
        U_plus = U*tmp.Vt';   #U*V_star;
        D_plus = Diagonal(tmp.S .^-1);
        D_plus[findall(x -> x == Inf, D_plus)] .= 0;

        # compute gain K
        K = U_plus*(D_plus*(D_plus*(U_plus'*(F'*(sqrtVinv'*sqrtVinv)))));

        # update estimate
        m = a + K*(Y[t:t]-F*a);

        store_a[:,t] = a;
        store_m[:,t+1] = m;
        store_D_plus[:,:,t+1] = D_plus;
        store_U_plus[:,:,t+1] = U_plus;
    end


    m = store_m[:,T+1];
    D_plus = store_D_plus[:,:,T+1];
    U_plus = store_U_plus[:,:,T+1];

    theta = m + U_plus*(D_plus*randn(q2))
    store_theta[:,T+1] = theta;

    tmp = svd(W);
    Dw = tmp.S .^0.5;
    Dw = max.(Dw, epss);
    sqrtWinv = Diagonal(1 ./ Dw) * tmp.Vt;
    
    # backward
    for t in collect(T:-1:1)

        a = store_a[:,t];
        m = store_m[:,t];
        D_plus = store_D_plus[:,:,t];
        U_plus = store_U_plus[:,:,t];

        D_plus_inv = D_plus ^-1;
        D_plus_inv[findall(x -> x == Inf, D_plus_inv)] .= 0;

        tmp = svd(Matrix([sqrtWinv*G*U_plus; D_plus_inv]));
        U_sq = U_plus*tmp.Vt';
        D_sq = 1 ./tmp.S;
        D_sq[findall(a-> a == Inf, D_sq)] .= 0;
        D_sq = Diagonal(D_sq);

        h = m + U_sq*D_sq*D_sq*U_sq'*G'*sqrtWinv'*sqrtWinv*(theta-a)
        theta = h + U_sq*(D_sq*randn(q2))

        store_theta[:,t] = theta
    end

    return store_theta'
end


# Example 1: Simulated Local Linear Trend Model with Seasonality
function ffbs_gibbs_sampler(y, nsim, init_psi, prior_shape, prior_rate)
    T = size(y, 1);
    N = nsim;
    Y = add_dim(y);

    # model matrices
    # y_t = F_t * theta_t + v_t   v_t ~ N(0, V_t)
    # theta_t = G_t * theta_t-1 + w_t   v_t ~ N(0, W_t)
    G, F, W, V = local_trend_seasonal12()

    q = size(W, 1)
    m0 = zeros(q, 1);
    C0 = 1e7I(q);
    theta = zeros(T+1, q)

    # prior hyperparameters
    # psi E(psi_x) = , Var(psi_y)=
    a_psiy = prior_shape[1];
    b_psiy = prior_rate[1];

    a_psi1 = prior_shape[2];
    b_psi1 = prior_rate[2];
    a_psi2 = prior_shape[3];
    b_psi2 = prior_rate[3];
    a_psi3 = prior_shape[4];
    b_psi3 = prior_rate[4];

    psiy = init_psi[1];
    psi1 = init_psi[2];
    psi2 = init_psi[3];
    psi3 = init_psi[4];

    new_a_psiy = a_psiy + T/2;
    new_a_psi1 = a_psi1 + T/2;
    new_a_psi2 = a_psi2 + T/2;
    new_a_psi3 = a_psi3 + T/2;

    store_psi_y = zeros(N);
    store_psi_1 = zeros(N);
    store_psi_2 = zeros(N);
    store_psi_3 = zeros(N);
    store_theta = zeros(T+1, q, N)

    for i in 1:N
        # FFBS
        theta = svd_forward_backward(Y, G, F, W, V, C0, m0);

        # draw phi_y
        ytheta = Y-theta[2:end,:]*F'
        SS_y = (ytheta'*ytheta)[1];
        new_b_psiy = b_psiy + 0.5 * SS_y;
        psiy = rand(Gamma(new_a_psiy, 1/new_b_psiy));
        V[1, 1] = 1/psiy;

        # SS_theta
        Δtheta = theta[2:end,:]-theta[1:end-1,:]*G';

        # draw psi_1
        SS_theta1 = Δtheta[:,1]'*Δtheta[:,1];
        new_b_psi1 = b_psi1 + 0.5 * SS_theta1[1];
        psi1 = rand(Gamma(new_a_psi1, 1/new_b_psi1));
        W[1, 1] = 1/psi1

        # draw psi_2
        SS_theta2 = Δtheta[:,2]'*Δtheta[:,2];
        new_b_psi2 = b_psi2 + 0.5 * SS_theta2[1];
        psi2 = rand(Gamma(new_a_psi2, 1/new_b_psi2));
        W[2, 2] = 1/psi2;

        # draw psi_3
        SS_theta3 = Δtheta[:,3]'*Δtheta[:,3];
        new_b_psi3 = b_psi3 + 0.5 * SS_theta3[1];
        psi3 = rand(Gamma(new_a_psi3, 1/new_b_psi3));
        W[3, 3] = 1/psi3

        if rem(i, N/10) == 0.0
            comp_perc = i/N*100
            println(string("completion: ", comp_perc, "%"))
        end

        store_psi_y[i] = psiy;
        store_psi_1[i] =  psi1;
        store_psi_2[i] =  psi2;
        store_psi_3[i] =  psi3;
        store_theta[:,:,i] = theta;
    end
    return store_psi_y, store_psi_1, store_psi_2, store_psi_3, store_theta
end
#
# y = simulate(250, local_trend_seasonal12, 12, 1, 5, 5e2, 1e5)
# plot(y)
#
# ret_hyper(5e2, 1e5)
#
# prior_shape = [1e-5, 2.5e-4, 2.5, 1e4];
# prior_rate = [1e-5, 5e-5, 0.005, 0.1];
# psi_init = [100, 100, 100, 100];
# nsim = 5000;
# Random.seed!(10);
#
# @time psi_y, psi_1, psi_2, psi_3, theta = gibbs_sampler_2(y, nsim, psi_init, prior_shape, prior_rate);
#
# # MCMC Diagnostics
# rho = acf(psi_y, 1000)
# scatter(collect(1:size(rho, 1)), rho)
#
# plot(cumsum(psi_y) ./ collect(1.0:nsim))
# plot(cumsum(psi_1) ./ collect(1.0:nsim))
# plot(cumsum(psi_2) ./ collect(1.0:nsim))
# plot(cumsum(psi_3) ./ collect(1.0:nsim))
#
# effective_sample_size(psi_y)
# effective_sample_size(psi_1)
# effective_sample_size(psi_2)
# effective_sample_size(psi_3)
#
# # psiy_hat
# mean(psi_y[50001:end])
# # psi1_hat
# mean(psi_1[50001:end])
# # psi2_hat
# mean(psi_2[5001:end])
# # psi3_hat
# mean(psi_3[50001:end])


# Example 2: Nile River
# Local Level with unknown Variance
function gibbs_sampler_2(y, nsim, nburn)

    T = size(y, 1);
    G = ones(1, 1);
    F = ones(1, 1);
    W = ones(1, 1);
    V = ones(1, 1);
    m0 = zeros(1, 1);
    C0 = 1e7 * ones(1, 1);
    Y = add_dim(y);

    a1 = 2;
    b1 = 0.0001;
    a2 = 2;
    b2 = 0.0001;

    new_a1 = a1 +T/2;
    new_a2 = a2 +T/2;

    N = nsim+nburn;

    store_phi1 = [];
    store_phi2 = [];
    store_theta = [];

    V = 1;
    W = 1;

    for i in 1:N
        # FFBS
        theta = FFBS(Y, G, F, [W], [V], C0, m0);
        theta = theta[end:-1:1];

        # draw phi_1
        ytheta = Y-theta
        new_b1 = b1 + 0.5 * (ytheta'*ytheta)[]
        phi1 = rand(Gamma(new_a1, 1/new_b1));
        V = phi1^-1

        # draw phi_2
        Δtheta = theta[2:T]-theta[1:T-1]
        new_b2 = b2 + 0.5 * (Δtheta'*Δtheta)[]
        phi2 = rand(Gamma(new_a2, 1/new_b2));
        W = phi2^-1

        if i > nburn
            push!(store_phi1, phi1);
            push!(store_phi2, phi2);
            append!(store_theta, theta);
        end
    end
    return store_phi1, store_phi2, store_theta
end
#
# data_raw = CSV.read("Nile.csv", header = 0);
# y = map(x->parse(Float64,x), data_raw[2:end, 2]);
# plot(y)
#
#
# @time store_phi1, store_phi2, store_theta = gibbs_sampler(y, 20, 10);
# @time store_phi1, store_phi2, store_theta = gibbs_sampler(y, 20000, 10000);
#
# plot(cumsum(store_phi1)./collect(1:20000))
# plot(cumsum(store_phi2)./collect(1:20000))
#
# phi1_hat = mean(store_phi1);
# phi2_hat = mean(store_phi2);
# theta = reshape(store_theta, 100, 20000);
# theta_hat = mean(theta, dims=2)
#
# plot(y)
# plot!(theta_hat)
#


# Example 3: AirPassengers
# Local Linear Trend + Seasonality with unknown Variance
function gibbs_sampler_3(y, nsim, init_psi, prior_shape, prior_rate)
    T = size(y, 1);
    N = nsim;
    Y = add_dim(y);

    # model matrices
    # y_t = F_t * theta_t + v_t   v_t ~ N(0, V_t)
    # theta_t = G_t * theta_t-1 + w_t   v_t ~ N(0, W_t)
    G, F, W, V = local_trend_seasonal12()

    q = size(W, 1)
    m0 = zeros(q, 1);
    C0 = 1e7I(q);
    theta = zeros(T+1, q)

    # prior hyperparameters
    # psi E(psi_x) = , Var(psi_y)=
    a_psiy = prior_shape[1];
    b_psiy = prior_rate[1];

    a_psi1 = prior_shape[2];
    b_psi1 = prior_rate[2];
    a_psi2 = prior_shape[3];
    b_psi2 = prior_rate[3];
    a_psi3 = prior_shape[4];
    b_psi3 = prior_rate[4];

    psiy = init_psi[1];
    psi1 = init_psi[2];
    psi2 = init_psi[3];
    psi3 = init_psi[4];

    new_a_psiy = a_psiy + T/2;
    new_a_psi1 = a_psi1 + T/2;
    new_a_psi2 = a_psi2 + T/2;
    new_a_psi3 = a_psi3 + T/2;

    store_psi_y = zeros(N);
    store_psi_1 = zeros(N);
    store_psi_2 = zeros(N);
    store_psi_3 = zeros(N);
    store_theta = zeros(T+1, q, N)

    for i in 1:N
        # FFBS
        theta = svd_forward_backward(Y, G, F, W, V, C0, m0);

        # draw phi_y
        ytheta = Y-theta[2:end,:]*F'
        SS_y = (ytheta'*ytheta)[1];
        new_b_psiy = b_psiy + 0.5 * SS_y;
        psiy = rand(Gamma(new_a_psiy, 1/new_b_psiy));
        V[1, 1] = 1/psiy;

        # SS_theta
        Δtheta = theta[2:end,:]-theta[1:end-1,:]*G';

        # draw psi_1
        SS_theta1 = Δtheta[:,1]'*Δtheta[:,1];
        new_b_psi1 = b_psi1 + 0.5 * SS_theta1[1];
        psi1 = rand(Gamma(new_a_psi1, 1/new_b_psi1));
        W[1, 1] = 1/psi1

        # draw psi_2
        SS_theta2 = Δtheta[:,2]'*Δtheta[:,2];
        new_b_psi2 = b_psi2 + 0.5 * SS_theta2[1];
        psi2 = rand(Gamma(new_a_psi2, 1/new_b_psi2));
        W[2, 2] = 1/psi2;

        # draw psi_3
        SS_theta3 = Δtheta[:,3]'*Δtheta[:,3];
        new_b_psi3 = b_psi3 + 0.5 * SS_theta3[1];
        psi3 = rand(Gamma(new_a_psi3, 1/new_b_psi3));
        W[3, 3] = 1/psi3

        store_psi_y[i] = psiy;
        store_psi_1[i] =  psi1;
        store_psi_2[i] =  psi2;
        store_psi_3[i] =  psi3;
        store_theta[:,:,i] = theta;
    end
    return store_psi_y, store_psi_1, store_psi_2, store_psi_3, store_theta
end

# Import data
#data_raw = CSV.read("./bayesian_inference/AirPassengers.csv", header = 0);
y=[112., 118, 132, 129, 121, 135, 148, 148, 136, 119, 104, 118, 115, 126, 141, 135, 125, 149, 170, 170, 158, 133, 114, 140, 145, 150, 178, 163, 172, 178,
 199, 199, 184, 162, 146, 166, 171, 180, 193, 181, 183, 218, 230, 242, 209, 191, 172, 194, 196, 196, 236, 235, 229, 243, 264, 272, 237, 211, 180, 201,
 204, 188, 235, 227, 234, 264, 302, 293, 259, 229, 203, 229, 242, 233, 267, 269, 270, 315, 364, 347, 312, 274, 237, 278, 284, 277, 317, 313, 318, 374,
 413, 405, 355, 306, 271, 306, 315, 301, 356, 348, 355, 422, 465, 467, 404, 347, 305, 336, 340, 318, 362, 348, 363, 435, 491, 505, 404, 359, 310, 337,
 360, 342, 406, 396, 420, 472, 548, 559, 463, 407, 362, 405, 417, 391, 419, 461, 472, 535, 622, 606, 508, 461, 390, 432]
#y = map(x->parse(Float64,x), data_raw[2:end, 2]);
y = broadcast(log, y);
plot(y)

ret_hyper(50, 1e5)

prior_shape = [1e-3, 2.5e-2, 2.5e5, 1e4];
prior_rate = [1e-4, 5e-4, 0.5, 0.1];
psi_init = [17, 63, 0.189, 0.71];
nsim = 100_000;
Random.seed!(10);

@time psi_y, psi_1, psi_2, psi_3, theta = gibbs_sampler_3(y, nsim, psi_init, prior_shape, prior_rate);

# MCMC Diagnostics
rho = acf(psi_y, 1000);
scatter(collect(1:size(rho, 1)), rho)

plot(cumsum(psi_y) ./ collect(1.0:nsim))
plot(cumsum(psi_1) ./ collect(1.0:nsim))
plot(cumsum(psi_2) ./ collect(1.0:nsim))
plot(cumsum(psi_3) ./ collect(1.0:nsim))

effective_sample_size(psi_y)
effective_sample_size(psi_1)
effective_sample_size(psi_2)
effective_sample_size(psi_3)

G, F, W, V = local_trend_seasonal12()

y_hat=mean(theta,dims=3)[:,:,1]*F';
l=y_hat-2*std(theta,dims=3)[:,:,1]*F';
u=y_hat+2*std(theta,dims=3)[:,:,1]*F';


x= [1:1:144;]
plot(x,y,label="Observed")
plot!(x,y_hat[2:end,1],label="Filtered", legend=:bottomright)
plot!(x, l[2:end],fillrange=u[2:end], alpha=.2, label="+-2 std")
savefig("FFBS_AirPassengers.png")

seas_mean=mean(theta,dims=3)[2:end,3,1];
seas_std=std(theta,dims=3)[2:end,3,1];
l=seas_mean-2*seas_std;
u=seas_mean+2*seas_std;
plot(x, seas_mean, label = "Seasonality")
plot!(x, l, fillrange = u, alpha = 0.2, label = "+-2*std")
savefig("FFBS_AirPassengers_seasonality.png")

# psiy_hat
mean(psi_y[5001:end])
# psi1_hat
mean(psi_1[5001:end])
# psi2_hat
mean(psi_2[5001:end])
# psi3_hat
mean(psi_3[5001:end])
