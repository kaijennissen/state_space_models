#------------------------------------------------------------------------------
# Julia Code for Kalman Filter
# Petris, G. & Petrone, S. & Campagnoli, P. (2009).
# Dynamic Linear Models with R.
#------------------------------------------------------------------------------
using LinearAlgebra
using CSV
using SparseArrays
using Plots
using Statistics
#using Distributions
#using Random


add_dim(x::Array) = reshape(x, (size(x)...,1))

# State Space Equation
# y_t = F_t * x_t + vega_t  vega_t ~ N(0, V_t)
# x_t = G_t * x_t-1 + omega_t  omega_t ~ N(0, W_t)

function kalman_filter(Y, G, F, W, V, C0, m0)

    theta = [];
    T = size(y, 1);
    m = m0;
    C = C0;

    for i in 1:T

        # pred state
        a = G*m
        R = G*C*G'+W

        # pred observation
        f = F*a
        Q = F*R*F'+V

        # filtered
        e = Y[i:i,]-f
        m = a + R*F'*inv(Q)*e
        C = R-R*F'*inv(Q)*F*R

        #push!(theta, m[:,1])
        append!(theta, m[:,1])
    end

    return theta
end


function kalman_filter_recursive(Y, G, F, W, V, C0, m0, t)
    theta = [];

    if t > 0
        # pred state
        m, C, theta = kalman_filter_recursive(Y, G, F, W, V, C0, m0, t-1)
        a = G*m
        R = G*C*G'+W

        # pred observation
        f = F*a
        Q = F*R*F'+V

        # filtered
        e = Y[t:t,]-f
        m = a + R*F'*inv(Q)*e
        C = R-R*F'*inv(Q)*F*R

        append!(theta, m[:,1])

    elseif t == 0
        m = m0
        C = C0
    end

    return m, C, theta
end



data_raw = CSV.read("AirPassengers.csv", header = 0);
y = map(x->parse(Float64,x), data_raw[2:end, 2]);
plot(y)
plot(broadcast(log, y))

# Local Level
# y_t = F * x_t + vega_t  vega_t ~ N(0, V_t)
# x_t = G * x_t-1 + omega_t  omega_t ~ N(0, W_t)
G = ones(1, 1);
F = ones(1, 1);
W = ones(1, 1);
V = ones(1, 1);
m0 = 100 * ones(1, 1);
C0 = 10.0 * ones(1, 1);
Y = add_dim(y);

theta1 = kalman_filter(Y, G, F, W, V, C0, m0);
m, C, theta2 = kalman_filter_recursive(Y, G, F, W, V, C0, m0, 144);

plot(theta1)
plot!(y)

# Local Linear Trend
# y_t = F * x_t + vega_t  vega_t ~ N(0, V_t)
# x_t = G * x_t-1 + omega_t  omega_t ~ N(0, W_t)
Y = add_dim(y);
F = [1 0];
V = ones(1, 1);
G = [1 1; 0 1];
W = [10 0; 0 5];
m0 = add_dim([100; 0]);
C0 = [10 0; 0 10];


theta1 = kalman_filter(Y, G, F, W, V, C0, m0);
m, C, theta2 = kalman_filter_recursive(Y, G, F, W, V, C0, m0, 144);

plot(theta1)
plot!(theta2)
plot!(y)

# Local Linear Trend + Seasonality
# y_t = F * x_t + vega_t  vega_t ~ N(0, V_t)
# x_t = G * x_t-1 + omega_t  omega_t ~ N(0, W_t)
y_log = broadcast(log, y);
plot(y_log)
Y = add_dim(y_log);

# initial values
X = [ones(144) 1:144]
beta_hat = X'X\X'y_log
y_hat = X*beta_hat
e = y_log - y_hat
error = reshape(e, 12, 12)
seas = mean(error, dims=2)

plot(y_log)
plot!(y_hat)

# Kalman Filter
F = [1 0 1 zeros(1, 10)];
V = 1.3*1e-4*ones(1, 1);
G = zeros(13, 13);
G[1, 1] = 1;
G[1, 2] = 1;
G[2, 2] = 1;
G[3, 3:end] = -1*ones(11);
G[4:13, 3:12] = 1.0I(10);
W = zeros(13, 13);
W[1, 1] = 7*1e-4;
W[2, 2] = 3.5*1e-13;
W[3, 3] = 6.4*1e-5;
m0 = [4.81; 0.01; seas[2:end]];
#m0 = mean(reshape(Y, 12,12), dims=2).-mean(reshape(Y, 12,12));
C0 = 1e7*I(13);
#C0[1, 1] = 0.01;
#C0[2, 2] = 0.01;
#C0[3, 3] = 0.5;

#theta1 = kalman_filter(Y, G, F, W, V, C0, m0);
m, C, theta2 = kalman_filter_recursive(Y, G, F, W, V, C0, m0, 144)
theta1 = reshape(theta2, 13, 144)'

plot(Y)
plot!(theta1[:,1])
plot!(theta1[:,1]+theta1[:,3])
