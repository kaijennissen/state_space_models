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

function SUR(X)
    r, c = size(X);
    idi = reshape((kron((1:r), ones(c, 1) )), r*c);
    idj = collect(1:r*c);
    X_SUR = sparse(idi, idj, reshape(X',r*c));
    return X_SUR
end

function gibbs_sampler(y::AbstractArray{Float64, 2}, nsim::Int64, burnin::Int64)

    y0 = y[1:2, :];
    Y = y[4:end, :];

    tt = size(y, 1)::Int64;
    nn = size(y, 2)::Int64;
    TT = size(Y, 1)::Int64;
    qq = nn*(nn+1)::Int64;
    TTqq = TT*qq::Int64;

    Y = reshape(Y',:,1);

    # priors #-----------------------------------------------------------------
    # Omega_11
    nu01 = nn + 3.::Float64
    S01 = 1.0I(nn)

    # Omega_22
    DD = 5. * sparse(1.0I,qq,qq)::SparseMatrixCSC{Float64,Int64};
    DD_inv = 1.0/5 * sparse(1.0I,qq,qq)::SparseMatrixCSC{Float64,Int64};
    nu02 = 6. * ones(qq);
    S02 = 0.01 * ones(qq);

    # initial values #---------------------------------------------------------
    new_nu01 = nu01 + TT;
    new_nu02 = broadcast(+ , nu02/2, (TT-1)/2)[:];

    # Omega11
    Omega11 = cov(y);
    Omega11_inv = inv(Omega11);

    # H
    H1 = sparse(1.0I, TTqq, TTqq)::SparseMatrixCSC{Float64,Int64};
    H2 = [[zeros(qq, TTqq-qq); sparse(1.0I, TTqq-qq, TTqq-qq)] zeros(TTqq, qq)]::SparseMatrixCSC{Float64,Int64};
    H = (H1 - H2)::SparseMatrixCSC{Float64,Int64};
    HT = sparse_transpose(H)::SparseMatrixCSC{Float64,Int64};

    # S
    Omega22 = 0.01 * sparse(1.0I, qq, qq)::SparseMatrixCSC{Float64,Int64};
    Omega22_inv = 10.0 * sparse(1.0I, qq, qq)::SparseMatrixCSC{Float64,Int64};
    S = blockdiag(DD, kron(sparse(1.0I,TT-1,TT-1), Omega22))::SparseMatrixCSC{Float64,Int64};
    S_inv = blockdiag(DD_inv, kron(sparse(1.0I,TT-1,TT-1), Omega22_inv))::SparseMatrixCSC{Float64,Int64};

    # G
    G = SUR([ones(TT*nn,1) kron(y[3:end-1, :], ones(nn))]);
    GT = sparse_transpose(G);

    # initialize for storeage
    store_beta = zeros(TTqq, nsim);
    store_Omega11 = zeros(nn, nn,  nsim);
    store_Omega22 = zeros(qq, nsim);

    for isim in 1:nsim+burnin

        #S_inv = blockdiag(DD_inv, kron(sparse(I,TT-1,TT-1), Omega22_inv));
        S_inv = kron(sparse(1.0I,TT,TT), Omega22_inv)::SparseMatrixCSC{Float64,Int64};
        S_inv[1:qq,1:qq] = DD_inv::SparseMatrixCSC{Float64,Int64};
        K = HT * S_inv * H::SparseMatrixCSC{Float64,Int64};

        GGL = tril(GT * kron(sparse(1.0I, TT, TT), Omega11_inv) * G)::SparseMatrixCSC{Float64,Int64};
        GT_Omega11_inv_G = (GGL + sparse_transpose(GGL) - sparse_diagonal(GGL))::SparseMatrixCSC{Float64,Int64};
        GT_Omega11_inv_Y = (GT * (kron(sparse(1.0I, TT, TT), Omega11_inv) * Y))::Array{Float64,2};

        P = K + GT_Omega11_inv_G::SparseMatrixCSC{Float64,Int64};

        C = cholesky(P, perm=1:TTqq);
        L = sparse(C.L)::SparseMatrixCSC{Float64,Int64};
        beta_hat = L'\(L\GT_Omega11_inv_Y)::Array{Float64,2};
        beta = (beta_hat + L' \ rand(Normal(),TTqq))::Array{Float64,2};

        # Omega11
        e1 = reshape(Y - G * beta, nn,:);
        new_S01 = S01 + e1*e1';
        Omega11 = rand(InverseWishart(new_nu01, new_S01));
        Omega11_inv = sparse(Symmetric(Omega11\I(nn)));

        # Omega22
        e2 = reshape(H * beta, qq, TT)';
        new_S02 = (S02 + sum(e2[2:end,:].^2, dims=1)[:])/2;
        for i in 1:qq
            Omega22[i,i] = rand(InverseGamma(new_nu02[i], new_S02[i]));
        end
        Omega22_inv = Omega22\sparse(1.0I,qq,qq);

        # store
        if isim > burnin
            i = isim - burnin;
            store_beta[:, i] = beta;
            store_Omega11[:,:, i] = Omega11;
            store_Omega22[:, i] = diag(Omega22);
        end
    end

    return store_beta, store_Omega11, store_Omega22
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
