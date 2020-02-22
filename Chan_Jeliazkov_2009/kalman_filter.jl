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
function SUR(X)
    r, c = size(X);
    idi = reshape((kron((1:r), ones(c, 1) )), r*c);
    idj = collect(1:r*c);
    X_SUR = sparse(idi, idj, reshape(X',r*c));
    return X_SUR
end

function gibbs_sampler(y::AbstractArray{Float64, 2}, nsim::Int64)

    epss = eps()^0.4

    y0 = y[1:2, :];
    Y = y[4:end, :];

    t = size(y, 1);
    n = size(y, 2);
    T = size(Y, 1);
    q = n*(n+1);
    Tq = T*q;

    Y = reshape(Y',:,1);

    # priors #-----------------------------------------------------------------
    # Omega_11
    nu01 = n + 3.;
    S01 = 1.0I(n);

    # Omega_22
    DD = 5. * sparse(1.0I, q, q);
    DD_inv = 1.0/5 * sparse(1.0I, q, q);
    nu02 = 6. * ones(q);
    S02 = 0.01 * ones(q);

    # initial values #---------------------------------------------------------
    new_nu01 = nu01 + T;
    new_nu02 = nu02/2 .+ (T-1)/2;

    # Omega11
    Omega11 = cov(y);
    Omega11_inv = inv(Omega11);

    # H
    H1 = sparse(1.0I, Tq, Tq);
    H2 = [[zeros(q, Tq-q); sparse(1.0I, Tq-q, Tq-q)] zeros(Tq, q)];
    H = (H1 - H2);
    Ht = sparse(H');

    # S
    Omega22 = 0.01 * sparse(1.0I, q, q);
    Omega22_inv = 10.0 * sparse(1.0I, q, q);
    S = blockdiag(DD, kron(sparse(1.0I, T-1, T-1), Omega22));
    S_inv = blockdiag(DD_inv, kron(sparse(1.0I, T-1, T-1), Omega22_inv));

    # G
    bigG = SUR([ones(T*n,1) kron(y[3:end-1, :], ones(n))]);
    bigGt = sparse(bigG');

    # initialize for storeage
    store_beta = zeros(Tq, nsim);
    store_Omega11 = zeros(n, n,  nsim);
    store_Omega22 = zeros(q, nsim);

    for isim in 1:nsim

        S_inv = blockdiag(DD_inv, kron(sparse(I,T-1,T-1), Omega22_inv));
        K = Ht * S_inv * H;

        GtOmega11inv = bigGt * kron(sparse(1.0I, T, T), Omega11_inv);
        GtOmega11invG = GtOmega11inv*bigG;
        GtOmega11invY = GtOmega11inv * Y;

        P = K + GtOmega11invG;
        P = (P+sparse(P'))./2

        tmp = cholesky(P, perm = 1:Tq, check = false);
        if issuccess(tmp)
            L = sparse(tmp.L);
        else
            tmp = ldlt(P, perm = 1:Tq, check = false);
            if issuccess(tmp)
                LD = sparse(tmp.LD);
                e = diag(LD);
                e[findall(x -> x <= epss, e)] .= 0;
                L = LD - Diagonal(LD) + I(Tq);
                D = sparse(Diagonal(e));
                L*D*sparse(L') ≈    P
                L = L*sparse(Diagonal(e .^0.5))
            end

        end

        beta_hat = sparse(L')\(L\GtOmega11invY);
        beta = beta_hat + sparse(L') \ rand(Normal(),Tq);

        # Omega11
        e1 = reshape(Y - bigG * beta, n,:);
        new_S01 = S01 + e1*e1';
        Omega11 = rand(InverseWishart(new_nu01, new_S01));
        Omega11_inv = sparse(Symmetric(Omega11\I(n)));

        # Omega22
        e2 = reshape(H * beta, q, T)';
        new_S02 = (S02 + sum(e2[2:end,:].^2, dims=1)[:])/2;
        for i in 1:q
            Omega22[i,i] = rand(InverseGamma(new_nu02[i], new_S02[i]));
        end
        Omega22_inv = Omega22\sparse(1.0I, q, q);

        # store
        store_beta[:, isim] = beta;
        store_Omega11[:,:, isim] = Omega11;
        store_Omega22[:, isim] = diag(Omega22);
    end

    return store_beta, store_Omega11, store_Omega22
end

# Kalman Filter
data_raw = CSV.read("./Chan_Jeliazkov_2009/USdata.csv", header = 0);
y = Matrix(data_raw);


@time store_beta, store_Omega11, store_Omega22  = gibbs_sampler(y, 10);
@time store_beta, store_Omega11, store_Omega22  = gibbs_sampler(y, Int(2.5e4));

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

savefig( "./Chan_Jeliazkov_2009/fg1")

# 2) Dynamic Factor Model
