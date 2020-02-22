# state space struct
using LinearAlgebra
using CSV
using SparseArrays
using Distributions
using Random
using Plots
using Documenter

# 1) State Space Model
# local level
# local linear trend
# seasonality

mutable struct StateSpaceModel
    # State Space Equation
    # y_t = X_t * beta_t + F_t * x_t + vega_t  vega_t ~ N(0, V_t)
    # eta_t = Z_t * gamma_t + G_t * eta_t-1 + omega_t  omega_t ~ N(0, W_t)

    y::AbstractArray{Float64, 2}
    t::Int
    n::Int
    q::Int
    T::Int
    Tq::Int

    F::AbstractArray{Float64, 2}

    G::AbstractArray{Float64, 2}

    X::AbstractArray{Float64, 2}

    Z::AbstractArray{Float64, 2}

    V::AbstractArray{Float64, 2}

    W::AbstractArray{Float64, 2}

    # constructor
    function StateSpaceModel(y, frequency = 0)
        T = size(y, 1);
        n = size(y, 2);
        q = np*2+n;
        Tq = T*q;
        bigY= reshape(Y',:,1);

        # Regular Model
        d = 2

        # Seasonality
        if frequency > 0
            ss = frequency - 1;

            # F
            Fs = [1 zeros(1, ss-1)];
            F = [F Fs]

            # W
            Ws = 1.0I(ss);
            W11 = W;
            W12 = zeros(d, ss);
            W21 = W12';
            W22 = Ws;
            W = [[W11 W12]; [W21 W22]]

            # G
            Gs = Matrix([-1.0*ones(1, ss);[1.0I(ss-1) zeros(ss-1)]]);
            G11 = G;
            G12 = zeros(d, ss);
            G21 = G12';
            G22 = Gs;
            G = [[G11 G12]; [G21 G22]]
        end

        new(y0, Y, bigY, t, n, q, T, Tq)
    end
end


mutable struct StateSpaceModel

    y0::AbstractArray{Float64, 2}
    Y::AbstractArray{Float64, 2}
    bigY::AbstractArray{Float64, 2}
    t::Int
    n::Int
    q::Int
    T::Int
    Tq::Int

    F

    G

    X

    Z

    V

    W



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




function construct(ssm)
    bigX = construct_bigX(ssm.TT, ssm.nn, ssm.qq, ssm.X);
    bigZ = construct_bigZ(ssm.TT, ssm.nn, ssm.qq, ssm.Z);
    bigGT = sparse_transpose(bigG);
    bigHT = sparse_transpose(bigH);
    bigG = construct_constant_bigG(ssm.TT, ssm.G);
    bigH = construct_constant_bigH(ssm.TT, ssm.qq, ssm.F);
    return bigX, bigZ, bigGT, bigHT, bigG, bigH
end










function sample_eta(ssm)
    bigS_inv = construct_constant_bigS_inv(DD_inv, Omega22_inv, TT);
    bigOmega11_inv = kron(sparse(1.0I,TT,TT), Omega11_inv);
    bigGamma = zeros(1,1);
    bigBeta = zeros(1, 1);
    K = bigHT * bigS_inv * bigH
    P = K + bigGT * bigOmega11_inv * bigG
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



function gibbs_sampler(ssm, nsim, burnin)

    bigX, bigZ, bigGT, bigHT, bigG, bigH = construct(ssm);

    # initialize storeage
    store_eta = zeros(ssm.TTqq, nsim);
    store_Omega11 = zeros(ssm.nn, ssm.nn,  nsim);
    store_Omega22 = zeros(ssm.qq, nsim);

    for i in 1:(burnin+nsim)

        # sample eta
        bigS_inv = construct_constant_bigS_inv(DD_inv, Omega22_inv, TT);
        bigOmega11_inv = kron(sparse(1.0I,TT,TT), Omega11_inv);
        bigGamma = zeros(1,1);
        bigBeta = zeros(1, 1);
        K = bigHT * bigS_inv * bigH
        P = K + bigGT * bigOmega11_inv * bigG
        C = cholesky(P, perm=1:ssm.TTqq);
        L = sparse(C.L);

        # eta_hat
        eta_tilde = bigH\(bigZ*bigGamma)
        eta_hat = K * eta_tilde + bigGT * (bigOmega11_inv * (ssm.bigY - ssm.bigX * ssm.bigBeta))

        # sample eta
        eta = (eta_hat + L' \ rand(Normal(),ssm.TTqq));

        # sample other variables
        # sample Omega11
        e1 = reshape(ssm.bigY - ssm.bigG * ssm.bigEta, ssm.nn,:);
        new_S01 = ssm.S01 + e1*e1';
        Omega11 = rand(InverseWishart(ssm.new_nu01, new_S01));
        Omega11_inv = sparse(Symmetric(Omega11\I(ssm.nn)));

        # sample Omega22
        e2 = reshape(ssm.bigH * ssm.bigEta, ssm.qq, ssm.TT)';
        new_S02 = (ssm.S02 + sum(e2[2:end,:].^2, dims=1)[:])/2;
        for i in 1:ssm.qq
            Omega22[i,i] = rand(InverseGamma(ssm.new_nu02[i], new_S02[i]));
        end
        Omega22_inv = Omega22\sparse(1.0I, ssm.qq, ssm.qq);

        # store
        if isim > burnin
            i = isim - burnin;
            store_eta[:, i] = eta;
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
