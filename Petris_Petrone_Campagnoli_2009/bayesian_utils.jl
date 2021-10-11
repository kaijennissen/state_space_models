# Collection on Utilities for Bayesian Analysis
using LinearAlgebra
#using CSV
#using Plots
#using Statistics
#using Distributions
#using Random
#using Infiltrator
#using JuliaInterpreter
#using BenchmarkTools

"""
Paramters
  ≡≡≡≡≡≡≡≡≡≡
n - dimension of measurement / observation vector
q - dimension of state vector
"""
function local_level()
    F = [1];
    G = [1];
    V = [1];
    W = [1];
    return G, F, W, V
end

add_dim(x::Array) = reshape(x, (size(x)...,1))

function local_trend_seasonal12()
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

function local_trend_seasonal4()
    # model
    F = [1.0 0.0 1.0 zeros(1, 2)];
    G = zeros(5, 5);
    G[1, 1] = 1.0;
    G[1, 2] = 1.0;
    G[2, 2] = 1.0;
    G[3, 3:end] = -1.0*ones(3);
    G[4:5, 3:4] = 1.0I(2);

    V = 1.0*ones(1, 1);
    W = zeros(5,5);
    W[1,1] = 1.0;
    W[2,2] = 1.0;
    W[3,3] = 1.0;

    return G, F, W, V
end

# ACF
function acf(y::Array{Float64,1}, k::Int64)
    rho_0 = var(y)
    acf = zeros(k)

    for i in 1:k
        acf[i] = cov(y[1:end-i],y[i+1:end])
    end

    return acf ./ rho_0
end

# Effective Sample Size
function effective_sample_size(x)
    rho = acf(x, size(x, 1)-2)
    n = size(x, 1)
    ess = n / (1+2*sum(rho))

    return ess
end

function ret_hyper(x, y)
    # alpha
    function ret_alpha(x, y)
         return x^2/y
    end

    # beta
    function ret_beta(x, y)
         return x/y
    end

    return ret_alpha(x, y), ret_beta(x, y)
end

# Simulate
function simulate(n, mod, freq, psiy, psi1, psi2, psi3)

    G, F, W, V = mod()

    V[1,1] = 1 ./psiy;
    W[1,1] = 1.0 / psi1;
    W[2,2] = 1.0 / psi2;
    W[3,3] = 1.0 / psi3;

    seas =  3*sin.(collect(range(-1, 1, length=freq)))
    theta0 =  [4 ; 0.1; seas[1:end-1]]

    # simulate
    y = zeros(n)

    theta =  theta0
    tmp = svd(W)
    sqrtW = tmp.U *Diagonal(tmp.S .^0.5)
    tmp = svd(V)
    sqrtV = tmp.U *(tmp.S .^0.5)

    for i in 1:n
        theta = G*theta+sqrtW*randn(size(sqrtW, 1))
        y[i] = (F*theta + sqrtV[]*randn(size(sqrtV, 1)))[]
    end
    return y
end
