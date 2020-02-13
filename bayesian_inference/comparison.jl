#------------------------------------------------------------------------------
# Comparison between traditional SVD-Forward-Filtering-Backward-Sampling and
# Kalman Filter of Chan & Jeliazkov 2009
#------------------------------------------------------------------------------
using LinearAlgebra
using CSV
using SparseArrays
using Distributions
using Random
using Plots
using Infiltrator

include("bayesian_inference/bayesian_utils.jl")
include("../Chan_Jeliazkov_2009/kalman_state_space_seasonal.jl")

y = simulate(500, local_trend_seasonal12, 12, 10, 50, 5e5, 1e5)
plot(y)

ret_hyper(50, 1e3)
# Chan & Jeliazkov 2009
nsim=10000
@time store_eta, store_Omega11_inv, store_Omega22_inv = gibbs_sampler(y, nsim)

# SVD-FFBS
prior_shape = [1e-3, 2.5e-2, 2.5e6, 2.5e6];
prior_rate = [1e-4, 5e-4, 5.0, 1.0];
psi_init = [100, 100, 100, 100];
@time psi_y, psi_1, psi_2, psi_3, theta = gibbs_sampler_1(y, nsim, psi_init, prior_shape, prior_rate);

plot(cumsum(store_Omega11_inv) ./ collect(1.0:nsim))
plot(cumsum(psi_y) ./ collect(1.0:nsim))
effective_sample_size(store_Omega11_inv)
effective_sample_size(psi_y)

rho = acf(store_Omega11_inv, 100)
scatter(collect(1:size(rho, 1)), rho)
rho = acf(psi_y, 100)
scatter(collect(1:size(rho, 1)), rho)

plot(cumsum(store_Omega22_inv[1,:]) ./ collect(1.0:nsim))
plot(cumsum(psi_1) ./ collect(1.0:nsim))
effective_sample_size(store_Omega22_inv[1,:])
effective_sample_size(psi_1)

rho = acf(store_Omega22_inv[1,:], 100)
scatter(collect(1:size(rho, 1)), rho)
rho = acf(psi_1, 100)
scatter(collect(1:size(rho, 1)), rho)
