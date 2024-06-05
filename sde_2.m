% Ronan Bottoms
% AMATH 536
% June 5, 2024


% Parameters for the model
% Currently set to Dan-G fit parameters
l = 0.139;
a = 0.473;
b = 0.206;
u = 15.54;
g = 0.069;
sat = 1000000000000;

% Adjustable growth rate parameters
ls = 0.139;
ld = 5.55*10^(-3);
r = 1;

duration = 100; % days
total_dosage = 1;
n0 = [10000000 0 0]; % 10^6 healthy cells
tspan = linspace(0, duration, 500);
global t0;
t0 = 0;

[t,n] = ode45(@(t,n) model(t,n,a,b,@I_133_uptake,total_dosage,u,g,l, ...
    ls,ld,sat, t0), tspan, n0);

% Plot simulation result
figure; hold on;
plot(t, n(:,1))
plot(t, n(:,2))
plot(t, n(:,3))


function delta_w = draw_incr(delta_t)
% Draws the next Wiener increments
    delta_w = zeros(1, 3);
    for i = 1:3
        delta_w(1) = normrnd(0, delta_t);
    end
end


function S = get_S(l, ls, ld, a, b, u, g, sat, r, n)
% Get the matrix S based on a vector of population sizes n.
    cov_n_n = l * n(1) * (1 - (n(1) + n(2) + n(3)) / sat) + u * n(2) +...
        (a+b)*r*n(1);
    cov_ns_ns = ls * n(2) * (1 - (n(1) + n(2) + n(3)) / sat) + u * n(2) + ...
        (a+b)*r*n(2) + b*r*n(1);
    cov_nd_nd = ld * n(3) * (1 - (n(1) + n(2) + n(3)) / sat) + (a+b)*r*n(2) +...
        g*n(3) + a*r*n(1);
    cov_ns_n = -u*n(2) - b*r*n(1);
    cov_nd_n = -a*r*n(1);
    cov_nd_ns = -(a+b)*r*n(2);
    
    V = [cov_n_n, cov_ns_n, cov_nd_n;
         cov_ns_n, cov_ns_ns, cov_nd_ns;
         cov_nd_n, cov_nd_ns, cov_nd_nd];
    
    S = sqrtm(V);
end


function dn = model(t, n, a, b, dosage, dose_amount, u, g, l, ls, ld, sat, t0)
% The full model. Dosage is a function that computes the amount of
% radiation in the system at time t.

r = dosage(t, dose_amount);
dn = zeros(3, 1);
% [N, Ns, Nd]
sat_term = (1 - (n(1) + n(2) + n(3)) / sat);
dn(1) = l * n(1) * sat_term + u * n(2) - (a + b) * r * n(1);
dn(2) = ls * n(2) * sat_term + b * r * n(1) - u * n(2) - (a + b) * r * n(2);
dn(3) = ld * n(3) * sat_term + a * r * (n(1) + n(2)) + b * r * n(2) - g * n(3);

% Compute stochastic element
S = get_S(l, ls, ld, a, b, u, g, sat, r, n);
dw = draw_incr(t - t0);
x = S * dw.';
dn = dn + x;

% Set previous time
t0 = t;
end


function dose_t = I_133_dose(t, amount)
% Returns an array of the dosage at time t (days) for a single injection
% of amount 'amount' of iodine 133

dose_t = amount * power(1/2, t / 8);
end


function dose_t = I_133_uptake(t, amount)
% Returns an array of the dosage at time t (days) for a single injection
% of amount 'amount' of iodine 133 that follows an uptake model.
dose_t = amount - amount * power(1/2, t / 8);
end


function dose_t = constant_dose_rate(t, amount)
% Returns an array of a constant dose rate at each time in t.

dose_t = zeros(length(t), 1) + amount;
end