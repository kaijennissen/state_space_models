#------------------------------------------------------------------------------
# Julia Code for a fast Kalman Filter based on
# Chan, J.C.C. and Jeliazkov, I. (2009). Efficient Simulation and
# Integrated Likelihood Estimation in State Space Models,
# International Journal of Mathematical Modelling and Numerical
# Optimisation, 1, 101-120.
#------------------------------------------------------------------------------
using LinearAlgebra
using CSV
using SparseArrays
using Distributions
using Random
using Plots
using Infiltrator

include("../bayesian_inference/bayesian_utils.jl")
# 1) TVP-VAR
function sparse_transpose(X::SparseMatrixCSC{Float64,Int64})
    i, j, v_ij = findnz(X);
    Xt = sparse(j, i, v_ij);
    return Xt
end

function sparse_diagonal(X::SparseMatrixCSC{Float64,Int64})
    i, j, v_ij = findnz(X);
    diag_elem = i.==j;
    new_i = i[diag_elem];
    new_v_ij = v_ij[diag_elem]
    Xt = sparse(new_i, new_i, new_v_ij);
    return Xt
end

# function SUR(X)
#     r, c = size(X);
#     idi = reshape((kron((1:r), ones(c, 1) )), r*c);
#     idj = collect(1:r*c);
#     X_SUR = sparse(idi, idj, reshape(X',r*c));
#     return X_SUR
# end

y = simulate(250, local_trend_seasonal12, 12, 10, 50, 5e5, 1e5)
plot(y)

ret_hyper(1e8, 1e3)

nsim=10000
@time store_eta, store_Omega11, store_Omega22 = gibbs_sampler(y, nsim)

plot(cumsum(store_Omega11) ./ collect(1.0:nsim))
plot(cumsum(store_Omega22[1,:]) ./ collect(1.0:nsim))
plot(cumsum(store_Omega22[2,:]) ./ collect(1.0:nsim))
plot(cumsum(store_Omega22[3,:]) ./ collect(1.0:nsim))
plot(cumsum(store_Omega22[4,:]) ./ collect(1.0:nsim))
plot(cumsum(store_Omega22[5,:]) ./ collect(1.0:nsim))
plot(cumsum(store_Omega22[6,:]) ./ collect(1.0:nsim))

eta_hat = reshape(mean(store_eta, dims=2), 13, :)'
plot(y)
plot!(eta_hat[:,1])

function gibbs_sampler(y, nsim, prior_shape, prior_rate)

    Y = add_dim(y);
    T = size(Y, 1);
    q = 13;
    Tq = T*q;

    # priors #-----------------------------------------------------------------

    # Omega_11
    a_y = prior_shape[1];
    b_y = prior_rate[1];

    # Omega_22
    DD_inv = sparse(1.0I, q, q);
    a_psi = prior_shape[2:end];
    b_psi = prior_rate[2:end];

    new_a_y = a_y + (T-1)/2;
    new_a_psi = a_psi .+ (T-1)/2;

    # initial values #---------------------------------------------------------

    # Omega11
    Omega11_inv = 1.0I(1);

    # G
    G = [1 0 1 zeros(1, 10)];

    # F
    F = zeros(q, q);
    F[1, 1] = 1;
    F[1, 2] = 1;
    F[2, 2] = 1;
    F[3, 3:end] = -1*ones(11);
    F[4:end, 3:end-1] = 1.0I(10);

    # H
    H1 = sparse(1.0I, Tq, Tq)::SparseMatrixCSC{Float64,Int64};
    H2 = [[zeros(q, Tq-q); kron(sparse(1.0I, T-1, T-1), F)] zeros(Tq, q)];
    H = (H1 - H2)::SparseMatrixCSC{Float64,Int64};
    HT = sparse(H')::SparseMatrixCSC{Float64,Int64};

    # S
    Omega22_inv = sparse(1.0I, q, q)::SparseMatrixCSC{Float64,Int64};
    S_inv = blockdiag(DD_inv, kron(sparse(1.0I, T-1, T-1), Omega22_inv))::SparseMatrixCSC{Float64,Int64};

    # G
    bigG = kron(sparse(1.0I, T, T), G);
    bigGT = sparse(bigG');

    # storeage
    store_eta = zeros(Tq, nsim);
    store_Omega11_inv = zeros(nsim);
    store_Omega22_inv = zeros(q, nsim);

    for isim in 1:nsim

        #S_inv = blockdiag(DD_inv, kron(sparse(I,TT-1,TT-1), Omega22_inv));
        S_inv = kron(sparse(1.0I, T, T), Omega22_inv);
        S_inv[1:q, 1:q] = DD_inv;
        K = (HT * S_inv * H);

        GGL = tril(bigGT * kron(sparse(1.0I, T, T), Omega11_inv) * bigG);
        GT_Omega11_inv_G = GGL + sparse(GGL') - Diagonal(GGL);
        GT_Omega11_inv_Y = bigGT * (kron(sparse(1.0I, T, T), Omega11_inv) * Y);

        P = K + GT_Omega11_inv_G;

        C = cholesky(P, perm=1:Tq);
        L = sparse(C.L);
        eta_hat = L'\(L\GT_Omega11_inv_Y);
        eta = (eta_hat + L' \ rand(Normal(), Tq));

        # Omega11
        e1 = Y - bigG * eta;
        new_b_y = b_y + (e1[2:end,:]'*e1[2:end,:])[]/2;
        Omega11_inv = rand(Gamma(new_a_y, 1/new_b_y));
        #Omega11 = Omega11_inv ^-1;

        # Omega22
        e2 = reshape(H * eta, q, T);
        new_b_psi = b_psi + diag(e2[:,2:end]*e2[:,2:end]')./2;
        for i in 1:q
            Omega22_inv[i, i] = rand(Gamma(new_a_psi[i], 1/new_b_psi[i]));
        end
        #Omega22 = Omega22_inv\sparse(1.0I, q, q);

        if rem(isim, nsim/10) == 0.0
            comp_perc = isim/nsim*100
            println(string("completion: ", comp_perc, "%"))
        end

        # store
        store_eta[:, isim] = eta;
        store_Omega11_inv[isim] = Omega11_inv;
        store_Omega22_inv[:, isim] = diag(Omega22_inv);
    end

    return store_eta, store_Omega11_inv, store_Omega22_inv
end


# Kalman Filter
data_raw = CSV.read("USdata.csv", header = 0);
y = Matrix(data_raw);

@time store_beta, store_Omega11, store_Omega22  = gibbs_sampler(y, 10, 10);
@time store_beta, store_Omega11, store_Omega22  = gibbs_sampler(y, 500, 500);
#@time store_beta, store_Omega11, store_Omega22  = gibbs_sampler(y, 50000, 5000);

beta_hat = mean(store_beta, dims = 2)';
Omega11_hat = mean(store_Omega11, dims = 3)[:,:,1];
Omega22_hat = mean(store_Omega22, dims = 2)';

beta = reshape(beta_hat, 4, 5, :);

l = @layout [a b ; c d];
p1 = plot(1:245, beta[1,:,:]', legend=false);
p2 = plot(1:245, beta[2,:,:]', legend=false);
p3 = plot(1:245, beta[3,:,:]', legend=false);
p4 = plot(1:245, beta[4,:,:]', legend=false);
plot(p1, p2, p3, p4, layout = l)

savefig( "fg1")

# 2) Dynamic Factor Model
