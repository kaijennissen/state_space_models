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

y = simulate(250, local_trend_seasonal12, 12, 1, 5, 5e2, 1e5)
plot(y)

ret_hyper(1e5, 1e5)

prior_shape = [1e-5, 2.5e-4, 2.5, 1e4, ];
prior_rate = [1e-5, 5e-5, 0.005, 0.1];
psi_init = [100, 100, 100, 100];
nsim = 5000;
Random.seed!(10);
nsim=1000
@time store_eta, store_Omega11, store_Omega22 = gibbs_sampler(y, nsim)

psi_y = store_Omega11
plot(cumsum(psi_y) ./ collect(1.0:nsim))
plot(cumsum(store_Omega22[2,:]) ./ collect(1.0:nsim))



function gibbs_sampler(y, nsim::Int64)

    Y = add_dim(y);
    tt = size(y, 1)::Int64;
    TT = size(Y, 1)::Int64;
    qq = 13;
    TTqq = TT*qq::Int64;


    # priors #-----------------------------------------------------------------
    # Omega_11
    a_y = 1e-5;
    b_y = 1e-5;

    # Omega_22
    DD = 5. * sparse(1.0I, qq, qq)::SparseMatrixCSC{Float64,Int64};
    DD_inv = 1.0/5 * sparse(1.0I,qq,qq)::SparseMatrixCSC{Float64,Int64};
    #psi_prior_shape = [2.5e-4, 2.5, 1e4ones(11,1)]
    #psi_prior_rate = [1e-5, 5e-5, 0.005, 1e2ones(11,1)];
    a_psi = [2.5e-4; 2.5; 1e5; 1e-17ones(10)];
    b_psi = [5e-5; 0.005; 1.0; 1e-9ones(10)];

    # initial values #---------------------------------------------------------
    new_a_y = a_y + TT/2;
    new_a_psi = a_psi .+ TT/2;

    # Omega11
    Omega11 = 1.0I(1);
    Omega11_inv = inv(Omega11);

    # G
    G = [1 0 1 zeros(1, 10)];

    # F
    F = zeros(qq, qq);
    F[1, 1] = 1;
    F[1, 2] = 1;
    F[2, 2] = 1;
    F[3, 3:end] = -1*ones(11);
    F[4:end, 3:end-1] = 1.0I(10);

    # H
    H1 = sparse(1.0I, TTqq, TTqq)::SparseMatrixCSC{Float64,Int64};
    H2 = [[zeros(qq, TTqq-qq); kron(sparse(1.0I, TT-1, TT-1), F)] zeros(TTqq, qq)];
    H = (H1 - H2)::SparseMatrixCSC{Float64,Int64};
    HT = sparse_transpose(H)::SparseMatrixCSC{Float64,Int64};

    # S
    #psi = ones(qq);
    #Omega22 = Diagonal(psi);
    Omega22 = sparse(1.0I, qq, qq)::SparseMatrixCSC{Float64,Int64};
    Omega22_inv = sparse(1.0I, qq, qq)::SparseMatrixCSC{Float64,Int64};
    S = blockdiag(DD, kron(sparse(1.0I,TT-1,TT-1), Omega22))::SparseMatrixCSC{Float64,Int64};
    S_inv = blockdiag(DD_inv, kron(sparse(1.0I,TT-1,TT-1), Omega22_inv))::SparseMatrixCSC{Float64,Int64};

    # G
    bigG = kron(sparse(1.0I, TT, TT), G);
    bigGT = sparse(bigG'); #sparse_transpose(bigG);

    # initialize for storeage
    store_eta = zeros(TTqq, nsim);
    store_Omega11 = zeros(nsim);
    store_Omega22 = zeros(qq, nsim);

    for isim in 1:nsim

        #S_inv = blockdiag(DD_inv, kron(sparse(I,TT-1,TT-1), Omega22_inv));
        S_inv = kron(sparse(1.0I,TT,TT), Omega22_inv)::SparseMatrixCSC{Float64,Int64};
        S_inv[1:qq,1:qq] = DD_inv::SparseMatrixCSC{Float64,Int64};
        K = HT * S_inv * H::SparseMatrixCSC{Float64,Int64};

        GGL = tril(bigGT * kron(sparse(1.0I, TT, TT), Omega11_inv) * bigG)::SparseMatrixCSC{Float64,Int64};
        GT_Omega11_inv_G = (GGL + sparse(GGL') - Diagonal(GGL))#''::SparseMatrixCSC{Float64,Int64};
        GT_Omega11_inv_Y = (bigGT * (kron(sparse(1.0I, TT, TT), Omega11_inv) * Y))::Array{Float64,2};

        P = K + GT_Omega11_inv_G::SparseMatrixCSC{Float64,Int64};

        C = cholesky(P, perm=1:TTqq);
        L = sparse(C.L)::SparseMatrixCSC{Float64,Int64};
        eta_hat = L'\(L\GT_Omega11_inv_Y)::Array{Float64,2};
        eta = (eta_hat + L' \ rand(Normal(),TTqq))::Array{Float64,2};

        # Omega11
        e1 = Y - bigG * eta;
        new_b_y = b_y + (e1'*e1)[]/2;
        #@infiltrate
        #Omega11 = rand(Gamma(new_nu01, 1/new_S01));
        Omega11_inv = rand(Gamma(new_a_y, 1/new_b_y));
        Omega11 = Omega11_inv ^-1;

        # Omega22
        e2 = reshape(H * eta, qq, TT);
        e2*e2'
        new_b_psi = b_psi + diag(e2*e2')./2;
        for i in 1:qq
            Omega22_inv[i,i] = rand(Gamma(new_a_psi[i], 1/new_b_psi[i]));
        end
        Omega22 = Omega22\sparse(1.0I,qq,qq);

        # store
        store_eta[:, isim] = eta;
        store_Omega11[isim] = Omega11_inv;
        store_Omega22[:, isim] = diag(Omega22_inv);
    end

    return store_eta, store_Omega11, store_Omega22
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
