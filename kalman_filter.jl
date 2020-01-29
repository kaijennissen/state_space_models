#------------------------------------------------------------------------------
# Julia Code for Kalman Filter
# Petris, G. & Petrone, S. & Campagnoli, P. (2009).
# Dynamic Linear Models with R.
#------------------------------------------------------------------------------
using LinearAlgebra
using CSV
using SparseArrays
using Distributions
using Random
using Plots

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

        push!(theta, m[1,])
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

        push!(theta, m[1,])

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
m0 = 100* ones(1, 1);
C0 = 10.0 * ones(1, 1);
Y = add_dim(y);

theta1 = kalman_filter(Y, G, F, W, V, C0, m0);
m, C, theta2 = kalman_filter_recursive(Y, G, F, W, V, C0, m0, 144);

plot(theta1)
plot!(y)

# Local Linear Trend
# y_t = F * x_t + vega_t  vega_t ~ N(0, V_t)
# x_t = G * x_t-1 + omega_t  omega_t ~ N(0, W_t)
y_log = broadcast(log, y);
plot(y_log)
Y = add_dim(y_log);
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

# Local Level + Deterministic Seasonality
# y_t = F * x_t + vega_t  vega_t ~ N(0, V_t)
# x_t = G * x_t-1 + omega_t  omega_t ~ N(0, W_t)
F = [1 0];
V = ones(1, 1);
G = [1 1; 0 1];
W = [10 0; 0 5];
m0 = add_dim([100; 0]);
C0 = [10 0; 0 10];
Y = add_dim(y);

theta1 = kalman_filter(Y, G, F, W, V, C0, m0);
m, C, theta2 = kalman_filter_recursive(Y, G, F, W, V, C0, m0, 144)

plot(theta1)
plot!(theta2)
plot!(y)

# Local Level + Stochastic Seasonality
# y_t = F * x_t + vega_t  vega_t ~ N(0, V_t)
# x_t = G * x_t-1 + omega_t  omega_t ~ N(0, W_t)
F = [1 0];
V = ones(1, 1);
G = [1 1; 0 1];
W = [10 0; 0 5];
m0 = add_dim([100; 0]);
C0 = [10 0; 0 10];
Y = add_dim(y);

theta1 = kalman_filter(Y, G, F, W, V, C0, m0);
m, C, theta2 = kalman_filter_recursive(Y, G, F, W, V, C0, m0, 144)

plot(theta1)
plot!(theta2)
plot!(y)
