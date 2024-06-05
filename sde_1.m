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

n = [1000000, 0, 0];
r = 1;

num_times = 10000;
tspan = linspace(0, 100, num_times);
pop = zeros(num_times, 3);
delta_t = 100 / num_times;
w = draw_incr(100 / num_times);

for i = 1:length(tspan)
    t = tspan(i);
    % increment the Wiener processes
    w = w + draw_incr(delta_t);
    % Compute S
    S = get_S(l, ls, ld, a, b, u, g, sat, r, n);
    % Compute SW
    x = (S * w.').';
    % Compute change in populations
    % Deterministic portion
    dn = zeros(1, 3);
    sat_term = (1 - (n(1) + n(2) + n(3)) / sat);
    dn(1) = l * n(1) * sat_term + u * n(2) - (a + b) * r * n(1);
    dn(2) = ls * n(2) * sat_term + b * r * n(1) - u * n(2) - (a + b) * r * n(2);
    dn(3) = ld * n(3) * sat_term + a * r * (n(1) + n(2)) + b * r * n(2) - g * n(3);
    % Add stochastic
    dn = dn + x;
    n = n + dn;
    % Record
    pop(i,:) = n;
end

% Plot simulation result


function delta_w = draw_incr(delta_t)
% Draws the next Wiener increments
    delta_w = zeros(1, 3);
    for i = 1:3
        delta_w(1) = normrnd(0, delta_t);
    end
    disp(delta_w)
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
         cov_nd_n, cov_nd_ns, cov_nd_nd]
    
    S = sqrtm(V);
end