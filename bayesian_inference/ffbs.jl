#------------------------------------------------------------------------------
# Julia Code for a Forward-Filtering Backward-Sampling Algorithm
# Petris, G. & Petrone, S. & Campagnoli, P. (2009).
# Dynamic Linear Models with R.
#------------------------------------------------------------------------------
using LinearAlgebra
using CSV
using Plots
using Statistics
using Distributions
using Random
using Infiltrator
using JuliaInterpreter

add_dim(x::Array) = reshape(x, (size(x)...,1))

function svd_forward_backward(Y, G, F, W, V, C0, m0)
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

    store_xm[:,1] = m;
    store_a[:,1] = a;
    store_m[:,1] = m;
    store_D_plus[:,:,1] = D_plus;
    store_U_plus[:,:,1] = U_plus;

    for t in 1:T

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


    m = store_m[:,T+1];
    D_plus = store_D_plus[:,:,T+1];
    U_plus = store_U_plus[:,:,T+1];

    #@infiltrate
    theta = m + U_plus*(D_plus*randn(q2))
    store_theta[:,T+1] = theta;

    #@infiltrate
    tmp = svd(W); #tmp <- La.svd(mod$W,nu=0)
    Dw = tmp.S .^0.5; #Dw <- sqrt(tmp$d)
    Dw = max.(Dw, epss); #Dw <- pmax(Dw, eps)
    #sqrtW = Diagonal(Dw); #Dw.inv <- 1/Dw
    sqrtWinv = Diagonal(1 ./ Dw) * tmp.Vt;
    #sqrtWinv <- Dw.inv * tmp$vt # t()%*%() = W^(-1)

    for t in collect(T:-1:1)

        a = store_a[:,t];
        m = store_m[:,t];
        D_plus = store_D_plus[:,:,t];
        U_plus = store_U_plus[:,:,t];
        #D = Diagonal(store_D[:,:,t]);
        #U = store_U[:,:,t];

        #@infiltrate
        D_plus_inv = D_plus ^-1; #D.inv <- 1/mod$D.C[i,];
        D_plus_inv[findall(x -> x == Inf, D_plus_inv)] .= 0; #D.inv[abs(D.inv)==Inf] <- 0
        #D_plus_inv = Diagonal(D_plus_inv);

        #La.svd(rbind(sqrtWinv %*% mod$GG %*% mod$U.C[[i]],diag(x=D.inv,nrow=length(D.inv))), nu=0)
        tmp = svd(Matrix([sqrtWinv*G*U_plus; D_plus_inv]));
        #@infiltrate
        U_sq = U_plus*tmp.Vt'; #U.H <- mod$U.C[[i]] %*% t(tmp$vt)
        D_sq = 1 ./tmp.S; #1/tmp$d; D.H[abs(D.H)==Inf] <- 0
        D_sq[findall(a-> a == Inf, D_sq)] .= 0;
        D_sq = Diagonal(D_sq);

        #@infiltrate
        h = m + U_sq*D_sq*D_sq*U_sq'*G'*sqrtWinv'*sqrtWinv*(theta-a)
        theta = h + U_sq*(D_sq*randn(q2))

        store_theta[:,t] = theta
    end

    #return store_m', store_xm'
    return store_theta'
end


# Example 1: Nile River
# Local Level with unknown Variance
data_raw = CSV.read("Nile.csv", header = 0);
y = map(x->parse(Float64,x), data_raw[2:end, 2]);
plot(y)

function gibbs_sampler_1(y, nsim, nburn)

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


@time store_phi1, store_phi2, store_theta = gibbs_sampler(y, 20, 10);
@time store_phi1, store_phi2, store_theta = gibbs_sampler(y, 20000, 10000);

plot(cumsum(store_phi1)./collect(1:20000))
plot(cumsum(store_phi2)./collect(1:20000))

phi1_hat = mean(store_phi1);
phi2_hat = mean(store_phi2);
theta = reshape(store_theta, 100, 20000);
theta_hat = mean(theta, dims=2)

plot(y)
plot!(theta_hat)



# Example 2: AirPassengers
# Local Linear Trend + Seasonality with unknown Variance
data_raw = CSV.read("./bayesian_inference/AirPassengers.csv", header = 0);
y = map(x->parse(Float64,x), data_raw[2:end, 2]);
y = broadcast(log, y);
plot(y);

function local_trend_seasonal()
    q = 13
    F = [1 0 1 zeros(1, 10)];
    G = zeros(q, q);
    G[1, 1] = 1;
    G[1, 2] = 1;
    G[2, 2] = 1;
    G[3, 3:end] = -1*ones(11);
    G[4:end, 3:end-1] = 1.0I(10);

    V = ones(1, 1);
    W = zeros(13,13);
    W[1,1] = 1;
    W[2,2] = 1;
    W[3,3] = 1;

    return G, F, W, V
end


function gibbs_sampler_2(y, nsim)
    T = size(y, 1);
    N = nsim;
    Y = add_dim(y);

    # model matrices
    # y_t = F_t * theta_t + v_t   v_t ~ N(0, V_t)
    # theta_t = G_t * theta_t-1 + w_t   v_t ~ N(0, W_t)
    G, F, W, V = local_trend_seasonal()

    q = size(W, 1)
    m0 = zeros(q, 1);
    C0 = 1e7I(q);
    # prior hyperparameters
    a_psiy = 2;
    b_psiy = 0.0001;

    a_psi1 = 2;
    b_psi1 = 0.0001;
    a_psi2 = 2;
    b_psi2 = 0.0001;
    a_psi3 = 2;
    b_psi3 = 0.0001;
    psiy = 1
    psi1 = 1
    psi2 = 1
    psi3 = 1

    new_a_psiy = a_psiy + T/2;
    new_a_psi1 = a_psi1 + T/2;
    new_a_psi2 = a_psi2 + T/2;
    new_a_psi3 = a_psi3 + T/2;

    store_psi_y = zeros(N);
    store_psi_1 = zeros(N);
    store_psi_2 = zeros(N);
    store_psi_3 = zeros(N);
    store_theta = zeros(T+1, 13, N)

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

@time gibbs_sampler_2(y, 100);

nsim = 30000
Random.seed!(10);
@time store_psi_y, store_psi_1, store_psi_2, store_psi_3, store_theta = gibbs_sampler_2(y, nsim);

plot(cumsum(store_psi_y) ./ collect(1.0:nsim))
plot(cumsum(store_psi_1) ./ collect(1.0:nsim))
plot(cumsum(store_psi_2) ./ collect(1.0:nsim))
plot(cumsum(store_psi_3) ./ collect(1.0:nsim))

psiy_hat = mean(store_psi_y[1001:end]);
psi1_hat = mean(store_psi_1[1001:end]);
psi2_hat = mean(store_psi_2[1001:end]);
psi3_hat = mean(store_psi_3[1001:end]);

theta_hat = mean(store_theta[:,:,5001:end], dims=3);
theta_max = maximum(store_theta[:,:,5001:end], dims=3);
theta_min = minimum(store_theta[:,:,5001:end], dims=3);
plot(y, legend=false)
plot!(theta_hat[2:end,1, 1])
plot!(theta_max[2:end,1, 1])
plot!(theta_min[2:end,1, 1])

# Diagnostics
tmp = zeros(144, nsim);
for i in 1:nsim
    tmp[:,i] = y-store_theta[2:end,:,i]*F';
end
SS_y[i] = tmp'*tmp;
theta_mc = reshape(theta_mc, 145,:)
theta_mc = broadcast(-,theta_mc[2:end,:],y)
SS_y = zeros(nsim);
for i in 1:nsim
    SS_y[i] = theta_mc[:,i]'*theta_mc[:,i]
end
plot(SS_y[101:end])


thetaG = zeros(13, 143, nsim)
for i in 1:nsim
  thetaG[:,:,i] = store_theta[:,2:end,i] - G*store_theta[:,1:end-1,i]
end

SS_theta1 = zeros(nsim)
for i in 1:nsim
  SS_theta1[i] = thetaG[1,:,i]'*thetaG[1,:,i]
end
plot(SS_theta1[101:end])
median(SS_theta1[101:end])



F = [1 0 1 zeros(1, 10)];
G = zeros(13, 13);
G[1, 1] = 1;
G[1, 2] = 1;
G[2, 2] = 1;
G[3, 3:end] = -1*ones(11);
G[4:13, 3:12] = 1.0I(10);
