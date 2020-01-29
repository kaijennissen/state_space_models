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
using Documenter

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

function SUR(X::AbstractArray{Float64, 2})
    r, c = size(X);
    idi = reshape((kron((1:r), ones(c, 1) )), r*c);
    idj = collect(1:r*c);
    X_SUR = sparse(idi, idj, reshape(X',r*c));
    return X_SUR
end

# bigX
function construct_bigX(TT, nn, qq, X)
    if isnothing(X)
        bigX = zeros(TT*nn, 1)
    else
        bigX = kron(ones(TT), X)
    end
    return bigX
end

# bigZ
function construct_bigZ(TT, nn, qq, Z)
    if isnothing(Z)
        bigZ = zeros(TT*qq, 1)
    else
        bigZ = kron(ones(TT), Z)
    end
    return bigZ
end

# bigG
function construct_constant_bigG(TT, G, dynamic = false)
    if dynamic
        bigG = blockdiag(G);
    else
        bigG = blockdiag(kron(sparse(1.0I, TT, TT), G));
    end
    return bigG
end

# bigH
function construct_constant_bigH(TT, qq, nn, F)
    TTqq = TT*qq;
    H1 = sparse(1.0I, TTqq, TTqq)::SparseMatrixCSC{Float64,Int64};
    H2 = [[zeros(qq, TTqq-qq); blockdiag(kron(sparse(1.0I, TT-1, TT-1), -F))] zeros(TTqq, qq)]::SparseMatrixCSC{Float64,Int64};
    bigH = (H1 + H2)::SparseMatrixCSC{Float64,Int64};
    return bigH;
end

# bigS
function construct_constant_bigS_inv(DD_inv, Omega22_inv, TT)
    S_inv = blockdiag(DD_inv, kron(sparse(1.0I,TT-1,TT-1), Omega22_inv));
    return S_inv;
end

# 1) State Space Model
mutable struct StateSpaceModel

    y0::AbstractArray{Float64, 2}
    Y::AbstractArray{Float64, 2}
    bigY::AbstractArray{Float64, 2}
    tt::Int
    nn::Int
    qq::Int
    TT::Int
    TTqq::Int

    # constructor
    function StateSpaceModel(y)
        y0 = y[1:3, :];
        Y = y[4:end, :];
        tt = size(y, 1);
        nn = size(y, 2);
        TT = size(Y, 1);
        qq = nn*(nn+1);
        TTqq = TT*qq;
        bigY= reshape(Y',:,1);
        new(y0, Y, bigY, tt, nn, qq, TT, TTqq)
    end
end



function sample_eta(ssm)
    ssm.bigS_inv = construct_constant_bigS_inv(DD_inv, Omega22_inv, TT);
    ssm.bigOmega11_inv = kron(sparse(1.0I,TT,TT), Omega11_inv);
    ssm.bigGamma = zeros(1,1);
    ssm.bigBeta = zeros(1, 1);
    K = ssm.bigHT * ssm.bigS_inv * ssm.bigH
    P = K + bigGT * ssm.bigOmega11_inv * ssm.bigG
    C = cholesky(P, perm=1:ssm.TTqq);
    L = sparse(C.L);

    # eta_hat
    eta_tilde = ssm.bigH\(ssm.bigZ*ssm.bigGamma)
    eta_hat = K * eta_tilde + ssm.bigGT * (ssm.bigOmega11_inv * (ssm.bigY - ssm.bigX * ssm.bigBeta))

    # sample eta
    eta = (eta_hat + L' \ rand(Normal(),ssm.TTqq));
    return eta
end

function sample_omega11(ssm)
    e1 = reshape(ssm.bigY - ssm.bigG * ssm.bigEta, ssm.nn,:);
    new_S01 = ssm.S01 + e1*e1';
    Omega11 = rand(InverseWishart(ssm.new_nu01, new_S01));
    Omega11_inv = sparse(Symmetric(Omega11\I(ssm.nn)));
    return Omega11
end

function sample_omega22(ssm)
    e2 = reshape(ssm.bigH * ssm.bigEta, ssm.qq, ssm.TT)';
    new_S02 = (ssm.S02 + sum(e2[2:end,:].^2, dims=1)[:])/2;
    for i in 1:ssm.qq
        Omega22[i,i] = rand(InverseGamma(ssm.new_nu02[i], new_S02[i]));
    end
    Omega22_inv = Omega22\sparse(1.0I, ssm.qq, ssm.qq);
    return Omega22_inv
end

function construct(ssm)
    ssm.bigX = construct_bigX(ssm.TT, ssm.nn, ssm.qq, ssm.X);
    ssm.bigZ = construct_bigZ(ssm.TT, ssm.nn, ssm.qq, ssm.Z);
    ssm.bigGT = sparse_transpose(ssm.bigG);
    ssm.bigHT = sparse_transpose(ssm.bigH);
    ssm.bigG = construct_constant_bigG(ssm.TT, ssm.G);
    ssm.bigH = construct_constant_bigH(ssm.TT, ssm.qq, ssm.F);
end


function gibbs_sampler(ssm, nsim, burnin)

    construct(ssm)

    # initialize storeage
    store_eta = zeros(ssm.TTqq, nsim);
    store_Omega11 = zeros(ssm.nn, ssm.nn,  nsim);
    store_Omega22 = zeros(ssm.qq, nsim);

    for i in 1:(burnin+nsim)

        # sample eta
        eta = sample_eta(ssm)

        # sample other variables
        # sample Omega11
        Omega11 = sample_omega11(ssm)

        # sample Omega22
        Omega2 = sample_omega22(ssm)

        # store
        if isim > burnin
            i = isim - burnin;
            store_beta[:, i] = beta;
            store_Omega11[:,:, i] = Omega11;
            store_Omega22[:, i] = diag(Omega22);
        end
    end

    return [store_beta, store_Omega11, store_Omega22]
end


# APPLICATION

data_raw = CSV.read("USdata.csv", header = 0);
y = Matrix(data_raw);

ssm = StateSpaceModel





















nn = 4;
qq = 20;
TT = 245;

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

# Omega11
Omega11 = cov(y);
Omega11_inv = inv(Omega11);

# Omega22
Omega22 = 0.01 * sparse(1.0I, qq, qq)::SparseMatrixCSC{Float64,Int64};
Omega22_inv = 10.0 * sparse(1.0I, qq, qq)::SparseMatrixCSC{Float64,Int64};


G = rand(nn, qq);
F = sparse(1.0I, qq, qq);







ssm = StateSpaceModel(y)






ssm = State_Space_Model(y, y, G, F, D, X = nothing, Z = nothing, DD, Omega11, Omega22)













function gibbs_sampler(y, nsim, burnin)
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
