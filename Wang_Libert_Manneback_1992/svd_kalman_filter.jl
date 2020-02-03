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
using JuliaInterpreter
using Infiltrator

function crossprod(X)
    return X'X
end

data_raw = CSV.read("./bayesian_inference/AirPassengers.csv", header = 0);
y = map(x->parse(Float64,x), data_raw[2:end, 2]);
plot(y)
Y=broadcast(log, y);

H = [1 0 1 zeros(1, 10)];
Psi = zeros(13, 13);
Psi[1, 1] = 1;
Psi[1, 2] = 1;
Psi[2, 2] = 1;
Psi[3, 3:end] = -1*ones(11);
Psi[4:13, 3:12] = 1.0I(10);

R = 1.0/22814.31*ones(1, 1);
Q = zeros(13, 13);
Q[1:3, 1:3] = 1.0*I(3);
Q[1, 1] = 1.0/2472.64;
Q[2, 2] = 1.0/55947.63;
Q[3, 3] = 1.0/4309.207;

G = 1.0I(13);

x0 = zeros(13,1);
P0 = 1e7I(13);

function svd_kalman(Y, H, Psi, G, R, Q, P0, x0)
    T = size(Y, 1);
    q1 = size(R, 1);
    q2 = size(Q, 1);

    store_a = zeros(q2, T);
    store_R = zeros(q2, q2, T);
    store_m = zeros(q2, T);
    store_C = zeros(q2, q2, T);

    store_x = zeros(q2, T+1)

    x̂ = x0;
    P = P0;

    UDUT = svd(P);
    U = Matrix(UDUT.U);
    D = Diagonal(UDUT.S);

    for t in 1:T

        # update U and D
        L = Matrix(cholesky(inv(R)).L)

        UDUT = svd(Matrix([L'*H*U ; D^-1]))
        V_star = UDUT.V
        D_star = Diagonal(UDUT.S)

        U_plus = U*V_star
        D_plus = D_star^-1

        # compute gain K
        K = U_plus*(D_plus^2)*U_plus'*H'*L*L'

        # update estimate x
        x̂_plus =  x̂ + K*(Y[t:t]-H*x̂)

        # project ahead
        Q_sq = Q^0.5
        UDUT = svd(Matrix([D_plus*U_plus'*Psi'; Q_sq'*G']))
        U = UDUT.V
        D = Diagonal(UDUT.S)

        x̂ = Psi*x̂_plus

        store_x[:,t] = x̂_plus;

    end
    return store_x
end

@time store_x = svd_kalman(Y, H, Psi, G, R, Q, C0, m0)

plot(1:145, store_x[1,:])

plot(Y)
plot!((H*store_x)[1:144])
plot!(([1 zeros(1,12)]*store_x)[1:144])
