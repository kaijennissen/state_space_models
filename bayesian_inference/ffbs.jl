#------------------------------------------------------------------------------
# Julia Code for a Forward-Filtering Backward-Sampling Algorithm
# Petris, G. & Petrone, S. & Campagnoli, P. (2009).
# Dynamic Linear Models with R.
#------------------------------------------------------------------------------
using LinearAlgebra
using CSV
using SparseArrays
using Plots
using Statistics


add_dim(x::Array) = reshape(x, (size(x)...,1))

# State Space Equation
# y_t = F_t * x_t + vega_t  vega_t ~ N(0, V_t)
# x_t = G_t * x_t-1 + omega_t  omega_t ~ N(0, W_t)

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

data_raw = CSV.read("AirPassengers.csv", header = 0);
y = map(x->parse(Float64,x), data_raw[2:end, 2]);
plot(y)
plot(broadcast(log, y))

# Local Level
# y_t = F * x_t + vega_t  vega_t ~ N(0, V_t)
# x_t = G * x_t-1 + omega_t  omega_t ~ N(0, W_t)
G = ones(1, 1);
F = ones(1, 1);
W = ones(1, 1);
V = ones(1, 1);
m0 = 100 * ones(1, 1);
C0 = 10.0 * ones(1, 1);
Y = add_dim(y);

theta = FFBS(Y, G, F, W, V, C0, m0)
Juno.@enter FFBS(Y, G, F, W, V, C0, m0);


function gibbs_sampler(y, nsim, nburn)

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


data_raw = CSV.read("Nile.csv", header = 0);
y = map(x->parse(Float64,x), data_raw[2:end, 2]);
plot(y)

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
