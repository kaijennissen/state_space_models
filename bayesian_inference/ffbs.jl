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

add_dim(x::Array) = reshape(x, (size(x)...,1))

function FFBS(Y, G, F, W, V, C0, m0)
    store_a = [];
    store_R = [];
    store_m = [];
    store_C = [];
    store_theta = []
    m = m0;
    C = C0;
    T = size(Y, 1);

    for t in 1:T
        # pred state
        a = G*m;
        R = G*C*G'+W;

        # pred observation
        f = F*a;
        Q = F*R*F'+V;

        # filtered
        e = Y[t:t,]-f;
        m = a + R*F'*inv(Q)*e;
        C = R-R*F'*inv(Q)*F*R;

        append!(store_a, a)
        append!(store_R, R)
        append!(store_m, m)
        append!(store_C, C)
    end

    h = store_m[T]
    H = store_C[T]
    #C = cholesky!(H);
    theta = h + sqrt(H)*randn()
    theta = add_dim([theta])
    append!(store_theta, theta)

    for t in collect((T-1):-1:1)
        a = add_dim([store_a[t+1]])
        R = add_dim([store_R[t+1]])
        m = add_dim([store_m[t]])
        C = add_dim([store_C[t]])

        h = m + C*G'*inv(R)*(theta-a)
        H = C - C*G'*inv(R)*G*C
        C = cholesky(H)
        theta = h + C.U*randn()
        append!(store_theta, theta)

    end
    return store_theta
end

function svd_forward_backward(Y, G, F, W, V, C0, m0)
    epss = eps(Float64)^0.4
    T = size(Y, 1);
    q1 = size(V, 1);
    q2 = size(W, 1);

    store_a = zeros(q2, T+1);
    store_m = zeros(q2, T+1);
    store_D_plus = zeros(q2, q2, T+1);
    store_U_plus = zeros(q2, q2, T+1);
    store_D = zeros(q2, q2, T);
    store_U = zeros(q2, q2, T);
    store_theta = zeros(q2, T+1)

    m = m0;
    C = C0;
    # predict
    a = G*m

    USVT = svd(C);
    U_plus = Matrix(USVT.U);
    D_plus = Matrix(Diagonal(USVT.S .^0.5));

    store_a[:,1] = a;
    store_m[:,1] = m;
    store_D_plus[:,:,1] = U_plus;
    store_U_plus[:,:,1] = U_plus;

    for t in 1:T



        # update U and D
        W_sq = W^0.5
        USVT = svd(Matrix([D_plus*U_plus'*G' ; W_sq]))
        U = USVT.V
        D = Diagonal(USVT.S)
        D_inv = Matrix(Diagonal(1 ./ USVT.S))

        L = Matrix(cholesky(inv(V)).L)

        USVT = svd(Matrix([L'*F*U; D_inv]))
        V_star = USVT.V
        D_star = Diagonal(USVT.S)
        D_star_inv = Diagonal(1 ./ USVT.S)

        U_plus = U*V_star
        D_plus = D_star_inv #D_star^-1

        # compute gain K
        K = U_plus*(D_plus^2)*U_plus'*F'*L*L'

        # update estimate x
        m = a + K*(Y[t:t]-F*a)
        # predict
        a = G*m

        store_a[:,t+1] = a;
        store_m[:,t+1] = m;
        store_D_plus[:,:,t+1] = D_plus;
        store_U_plus[:,:,t+1] = U_plus;
        #store_D[:,:,t] = D;
        #store_U[:,:,t] = U;

    end

    m = store_m[:,T+1];
    D_plus = store_D_plus[:,:,T+1];
    U_plus = store_U_plus[:,:,T+1];

    h = m;

    theta = h + U_plus*D_plus*randn(q2)
    store_theta[:,T+1] = theta;

    # W is time invariant
    tmp = svd(W)
    DW = tmp.S .^ 0.5
    DW = maximum([DW epss*ones(13,1)], dims=2)
    DW_inv = Diagonal(1 ./ DW)

    sqrtWinv = tmp.Vt'*DW_inv^2*tmp.Vt

    L_star = DW_inv*tmp.Vt

    for t in collect(T:-1:1)

        a = store_a[:,t+1];
        m = store_m[:,t];
        D_plus = store_D_plus[:,:,t];
        U_plus = store_U_plus[:,:,t];
        #D = store_D[:,:,t];
        #U = store_U[:,:,t];

        tmp = svd(Matrix([L_star'*G*U_plus; D_plus^-1]))
        V_tri = tmp.Vt

        U_sq = U_plus*V_tri'
        D_sq = tmp.S .^ -1

        D_sq[findall(a -> a == Inf, D_sq)] .= 0
        D_sq = Diagonal(D_sq)

        h = m + U_sq*D_sq^2*U_sq'*G'sqrtWinv*(theta-a)
        theta = h + U_sq*D_sq*randn(q2)

        store_theta[:,t] = theta

    end

    return store_theta
end

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
#plot(broadcast(log, y))

function gibbs_sampler_2(y, nsim)
    T = size(y, 1);
    N = nsim;
    Y = add_dim(y);

    # model matrices
    # y_t = F_t * theta_t + v_t   v_t ~ N(0, V_t)
    # theta_t = G_t * theta_t-1 + w_t   v_t ~ N(0, W_t)
    F = [1 0 1 zeros(1, 10)];
    G = zeros(13, 13);
    G[1, 1] = 1;
    G[1, 2] = 1;
    G[2, 2] = 1;
    G[3, 3:end] = -1*ones(11);
    G[4:13, 3:12] = 1.0I(10);

    V = ones(1, 1);
    W = zeros(13, 13);
    W[1:3, 1:3] = 1.0*I(3);

    m0 = zeros(13, 1);
    C0 = 1e7*I(13);

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
    store_theta = zeros(13, T+1, N)

    for i in 1:N
        # FFBS
        theta = svd_forward_backward(Y, G, F, W, V, C0, m0);
        #@infiltrate
        #theta = theta[:,end:-1:1];

        # draw phi_y
        ytheta = Y-(F*theta[:,2:end])'
        SS_y = (ytheta'*ytheta)[1];
        new_b_psiy = b_psiy + 0.5 * SS_y;
        psiy = rand(Gamma(new_a_psiy, 1/new_b_psiy));
        V[1, 1] = 1/psiy;

        # SS_theta
        Δtheta = theta[:,2:end]-G'theta[:,1:end-1];

        # draw psi_1
        SS_theta1 = Δtheta[1,:]*Δtheta[1,:]';
        new_b_psi1 = b_psi1 + 0.5 * SS_theta1[1];
        psi1 = rand(Gamma(new_a_psi1, 1/new_b_psi1));
        W[1, 1] = 1/psi1

        # draw psi_2
        SS_theta2 = Δtheta[2,:]*Δtheta[2,:]';
        new_b_psi2 = b_psi2 + 0.5 * SS_theta2[1];
        psi2 = rand(Gamma(new_a_psi2, 1/new_b_psi2));
        W[2, 2] = 1/psi2;

        # draw psi_3
        SS_theta3 = Δtheta[3,:]*Δtheta[3,:]';
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

@time gibbs_sampler_2(y, 3);

nsim = 1000
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

theta_hat = mean(store_theta, dims=3)[:,:,1];

plot(y, legend=false)
plot!(theta_hat[1,:])

# Diagnostics
theta_mc = reshape(store_theta, 13, :)'*F';
theta_mc = reshape(theta_mc, 145,:)
SS_y = zeros(nsim);
for i in 1:nsim
    SS_y[i] = theta_mc[:,i]'*theta_mc[:,i]
end
plot(SS_y[101:end])
median(SS_y[101:end])


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
