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

include("./bayesian_utils.jl")
include("../Chan_Jeliazkov_2009/kalman_state_space_seasonal.jl")
include("./ffbs.jl")

# compare using simulate data
y = simulate(1000, local_trend_seasonal12, 12, 10, 50, 5e5, 1e5)
plot(y)

# Chan & Jeliazkov 2009 -------------------------------------------------------
prior_shape = [1e-3; 2.5e-2; 2.5e8; 1.e7; 1e13ones(10)];
prior_rate = [1e-4; 5e-4; 500; 100; 1e5ones(10)];
nsim = 100000

@time eta, Omega11_inv, Omega22_inv = gibbs_sampler(y, nsim, prior_shape, prior_rate)

# Omega11_inv
x =  Omega11_inv;
nburn = 2500;
plot(cumsum(x) ./ collect(1.0:length(x)))
histogram(x[nburn:end])
effective_sample_size(x)
rho = acf(x, 100);
scatter(collect(1:size(rho, 1)), rho)

# Omega22_inv[1, 1]
x = Omega22_inv[1, :]
nburn = 100
plot(cumsum(x) ./ collect(1.0:length(x)))
histogram(x[nburn:end])
effective_sample_size(x)
rho = acf(x, 100)
scatter(collect(1:size(rho, 1)), rho)

#  x = Omega22_inv[2, 2]
x = Omega22_inv[2,:]
nburn = 100
plot(cumsum(x) ./ collect(1.0:length(x)))
histogram(x[nburn:end])
effective_sample_size(x)
rho = acf(x, 100)
scatter(collect(1:size(rho, 1)), rho)

#  Omega22_inv[3, 3]
x = Omega22_inv[3,:];
nburn = 100
plot(cumsum(x) ./ collect(1.0:length(x)))
histogram(x[nburn:end])
effective_sample_size(x);
rho = acf(x, 100);
scatter(collect(1:size(rho, 1)), rho);

# Omega2_inv[4, 4]
x = Omega22_inv[4,:]
nburn = 100
plot(cumsum(x) ./ collect(1.0:length(x)))
histogram(x[1001:end])
effective_sample_size(x)
rho = acf(x, 100)
scatter(collect(1:size(rho, 1)), rho)

# SVD-FFBS --------------------------------------------------------------------
prior_shape = [1e-3, 2.5e-2, 2.5e6, 2.5e6];
prior_rate = [1e-4, 5e-4, 5.0, 1.0];
psi_init = [100, 100, 100, 100];
nsim=10000
@time psi_y, psi_1, psi_2, psi_3, theta = ffbs_gibbs_sampler(y, nsim, psi_init, prior_shape, prior_rate);


# Diagnostics

# psi_y
x = psi_y
plot(cumsum(x) ./ collect(1.0:length(x)))
histogram(x[1001:end])
effective_sample_size(x)
rho = acf(x, 100)
scatter(collect(1:size(rho, 1)), rho)

# psi_1
x = psi_1
plot(cumsum(x) ./ collect(1.0:length(x)))
histogram(x[1001:end])
effective_sample_size(x)
rho = acf(x, 100)
scatter(collect(1:size(rho, 1)), rho)

# psi_2
x = psi_2
plot(cumsum(x) ./ collect(1.0:length(x)))
histogram(x[1001:end])
effective_sample_size(x)
rho = acf(x, 100)
scatter(collect(1:size(rho, 1)), rho)

# psi_3
x = psi_3
plot(cumsum(x) ./ collect(1.0:length(x)))
histogram(x[1001:end])
effective_sample_size(x)
rho = acf(x, 100)
scatter(collect(1:size(rho, 1)), rho)
