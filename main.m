%
% Ronan Bottoms
% AMATH 536
% June 5, 2024
%
% Contains functions for modeling the three types of models
% and creates plots for several different r values under multiple
% cell types and parameter regimes.


%plot_regen_rates()
plot_patu_dose()
plot_dan_g_dose()


function plot_patu_dose()
    % Radiation doses to enumerate over
    r = [10, 12, 15];
    
    figure;
    set(gcf,'units','normalized','outerposition',[0 0 0.3 0.7]);
    for i = 1:length(r)
        % General
        subplot(3, 2, i + (i-1)); hold on;
        [t, n] = Patu_8988T_gen_dose(0.139, 5.55*10^(-3), r(i), 100);
        plot(t, n(:, 1), "LineWidth", 2);
        plot(t, n(:, 2), "LineWidth", 2,'LineStyle','-.');
        plot(t, n(:, 3), "LineWidth", 2,'LineStyle',':');
        subtitle("r = " + r(i))
        ylabel("num cells")
    
        if (i == 1)
            title("Full Model")
        end
    
        if (i == 3)
            xlabel("time (days)")
        end
            
    
        % Simplified
        subplot(3, 2, 2 * i); hold on;
        [t, n] = Patu_8988T_simp_dose(r(i), 100);
        plot(t, n(:, 1), "LineWidth", 2);
        plot(t, n(:, 2), "LineWidth", 2,'LineStyle','-.');
        plot(t, n(:, 3), "LineWidth", 2,'LineStyle',':');
        subtitle("r = " + r(i))
        ylabel("num cells")
    
        if (i == 1)
            title("Simplified Model")
        end
    
        if (i == 3)
            xlabel("time (days)")
        end
    end
end

% Plot Dan-G Uptake
function plot_dan_g_dose()
% Radiation doses to enumerate over
r = [2, 2.9, 4];

figure;
set(gcf,'units','normalized','outerposition',[0 0 0.3 0.7]);
for i = 1:length(r)
    % General
    subplot(3, 2, i + (i-1)); hold on;
    [t, n] = Dan_G_gen_dose(0.139, 5.55*10^(-3), r(i), 100);
    plot(t, n(:, 1), "LineWidth", 2);
    plot(t, n(:, 2), "LineWidth", 2,'LineStyle','-.');
    plot(t, n(:, 3), "LineWidth", 2,'LineStyle',':');
    subtitle("r = " + r(i))
    ylabel("num cells")

    if (i == 1)
        title("Full Model")
    end

    if (i == 3)
        xlabel("time (days)")
    end
        

    % Simplified
    subplot(3, 2, 2 * i); hold on;
    [t, n] = Dan_G_simp_dose(r(i), 100);
    plot(t, n(:, 1), "LineWidth", 2);
    plot(t, n(:, 2), "LineWidth", 2,'LineStyle','-.');
    plot(t, n(:, 3), "LineWidth", 2,'LineStyle',':');
    subtitle("r = " + r(i))
    ylabel("num cells")

    if (i == 1)
        title("Simplified Model")
    end

    if (i == 3)
        xlabel("time (days)")
    end
end
end


% Plot Patu uptake
function plot_patu_uptake()
% Radiation doses to enumerate over
r = [1, 2, 3];

figure;
set(gcf,'units','normalized','outerposition',[0 0 0.3 0.7]);
for i = 1:length(r)
    % General
    subplot(3, 2, i + (i-1)); hold on;
    [t, n] = Patu_8988T_gen_uptake(0.139, 5.55*10^(-3), r(i), 100);
    plot(t, n(:, 1), "LineWidth", 2);
    plot(t, n(:, 2), "LineWidth", 2,'LineStyle','-.');
    plot(t, n(:, 3), "LineWidth", 2,'LineStyle',':');
    subtitle("r = " + r(i))
    ylabel("num cells")

    if (i == 1)
        title("Full Model")
    end

    if (i == 3)
        xlabel("time (days)")
    end
        

    % Simplified
    subplot(3, 2, 2 * i); hold on;
    [t, n] = Patu_8988T_simp_uptake(r(i), 100);
    plot(t, n(:, 1), "LineWidth", 2);
    plot(t, n(:, 2), "LineWidth", 2,'LineStyle','-.');
    plot(t, n(:, 3), "LineWidth", 2,'LineStyle',':');
    subtitle("r = " + r(i))
    ylabel("num cells")

    if (i == 1)
        title("Simplified Model")
    end

    if (i == 3)
        xlabel("time (days)")
    end
end
end

% Plot Dan-G Uptake
function plot_dan_g_uptake()
% Radiation doses to enumerate over
r = [0.1, 0.2, 0.4];

figure;
set(gcf,'units','normalized','outerposition',[0 0 0.3 0.7]);
for i = 1:length(r)
    % General
    subplot(3, 2, i + (i-1)); hold on;
    [t, n] = Dan_G_gen_uptake(0.139, 5.55*10^(-3), r(i), 100);
    plot(t, n(:, 1), "LineWidth", 2);
    plot(t, n(:, 2), "LineWidth", 2,'LineStyle','-.');
    plot(t, n(:, 3), "LineWidth", 2,'LineStyle',':');
    subtitle("r = " + r(i))
    ylabel("num cells")

    if (i == 1)
        title("Full Model")
    end

    if (i == 3)
        xlabel("time (days)")
    end
        

    % Simplified
    subplot(3, 2, 2 * i); hold on;
    [t, n] = Dan_G_simp_uptake(r(i), 100);
    plot(t, n(:, 1), "LineWidth", 2);
    plot(t, n(:, 2), "LineWidth", 2,'LineStyle','-.');
    plot(t, n(:, 3), "LineWidth", 2,'LineStyle',':');
    subtitle("r = " + r(i))
    ylabel("num cells")

    if (i == 1)
        title("Simplified Model")
    end

    if (i == 3)
        xlabel("time (days)")
    end
end
end

%% Plot regeneration rates
function plot_regen_rates()
a = 0.0;
b = 0.664;
lambda = 0.139;
gamma = 0.069;
sat = Inf;
n0 = [10000000 0 0]; % 10^6 healthy cells
tspan = linspace(0, 100, 500);
total_dosage = 2;

figure;
set(gcf,'units','normalized','outerposition',[0 0 0.5 0.2]);

subplot(1, 4, 1); hold on;
mu = 15.84;
[t,n] = ode45(@(t,n) model(t,n,a,b,@I_133_uptake,total_dosage,mu,gamma,lambda, ...
    0,0,sat), tspan, n0);

plot(t, n(:, 1), "LineWidth", 2);
plot(t, n(:, 2), "LineWidth", 2,'LineStyle','-.');
plot(t, n(:, 3), "LineWidth", 2,'LineStyle',':');
xlabel("time (days)");
subtitle("\mu = " + mu);
ylabel("num cells");

subplot(1, 4, 2); hold on;
mu = 7;
[t,n] = ode45(@(t,n) model(t,n,a,b,@I_133_uptake,total_dosage,mu,gamma,lambda, ...
    0,0,sat), tspan, n0);

plot(t, n(:, 1), "LineWidth", 2);
plot(t, n(:, 2), "LineWidth", 2,'LineStyle','-.');
plot(t, n(:, 3), "LineWidth", 2,'LineStyle',':');
xlabel("time (days)");
subtitle("\mu = " + mu);

subplot(1, 4, 3); hold on;
mu = 2;
[t,n] = ode45(@(t,n) model(t,n,a,b,@I_133_uptake,total_dosage,mu,gamma,lambda, ...
    0,0,sat), tspan, n0);

plot(t, n(:, 1), "LineWidth", 2);
plot(t, n(:, 2), "LineWidth", 2,'LineStyle','-.');
plot(t, n(:, 3), "LineWidth", 2,'LineStyle',':');
xlabel("time (days)");
subtitle("\mu = " + mu);

subplot(1, 4, 4); hold on;
mu = 1;
[t,n] = ode45(@(t,n) model(t,n,a,b,@I_133_uptake,total_dosage,mu,gamma,lambda, ...
    0,0,sat), tspan, n0);

plot(t, n(:, 1), "LineWidth", 2);
plot(t, n(:, 2), "LineWidth", 2,'LineStyle','-.');
plot(t, n(:, 3), "LineWidth", 2,'LineStyle',':');
xlabel("time (days)");
ylim([0 Inf]);
subtitle("\mu = " + mu);
end
%% Patu


function [t, n] = Patu_8988T_gen_uptake(lambda_s, lambda_d, total_dosage, duration)
a = 0.0;
b = 0.664;
lambda = 0.139;
gamma = 0.069;
mu = 0.011 * 1440; % convert from /min to /day to match lambdas
sat = Inf;
n0 = [10000000 0 0]; % 10^6 healthy cells
tspan = linspace(0, duration, 500);

[t,n] = ode45(@(t,n) model(t,n,a,b,@I_133_uptake,total_dosage,mu,gamma,lambda, ...
    lambda_s,lambda_d,sat), tspan, n0);
end

function [t, n] = Patu_8988T_simp_uptake(total_dosage, duration)
a = 0.0;
b = 0.664;
lambda = 0.139;
gamma = 0.069;
mu = 0.011 * 1440; % convert from /min to /day to match lambdas
sat = Inf;
n0 = [10000000 0 0]; % 10^6 healthy cells
tspan = linspace(0, duration, 500);

[t,n] = ode45(@(t,n) simplified_model(t,n,a,b,@I_133_uptake,total_dosage,mu,gamma,lambda, ...
    0,0,sat), tspan, n0);
end

function [t, n] = Patu_8988T_gen_dose(lambda_s, lambda_d, total_dosage, duration)
a = 0.0;
b = 0.664;
lambda = 0.139;
gamma = 0.069;
mu = 0.011 * 1440; % convert from /min to /day to match lambdas
sat = Inf;
n0 = [10000000 0 0]; % 10^6 healthy cells
tspan = linspace(0, duration, 500);

[t,n] = ode45(@(t,n) model(t,n,a,b,@I_133_dose,total_dosage,mu,gamma,lambda, ...
    lambda_s,lambda_d,sat), tspan, n0);
end

function [t, n] = Patu_8988T_simp_dose(total_dosage, duration)
a = 0.0;
b = 0.664;
lambda = 0.139;
gamma = 0.069;
mu = 0.011 * 1440; % convert from /min to /day to match lambdas
sat = Inf;
n0 = [10000000 0 0]; % 10^6 healthy cells
tspan = linspace(0, duration, 500);

[t,n] = ode45(@(t,n) simplified_model(t,n,a,b,@I_133_dose,total_dosage,mu,gamma,lambda, ...
    0,0,sat), tspan, n0);
end

%% Dan-G

function [t,n] = Dan_G_gen_dose(lambda_s, lambda_d, total_dosage, duration)
a = 0.473;
b = 0.206;
lambda = 0.139;
gamma = 0.069;
mu = 0.011 * 1440; % convert from /min to /day to match lambdas
sat = Inf;
n0 = [10000000 0 0]; % 10^6 healthy cells
tspan = linspace(0, duration, 500);

[t,n] = ode45(@(t,n) model(t,n,a,b,@I_133_dose,total_dosage,mu,gamma,lambda, ...
    lambda_s,lambda_d,sat), tspan, n0);
end

function [t,n] = Dan_G_simp_dose(total_dosage, duration)
a = 0.473;
b = 0.206;
lambda = 0.139;
gamma = 0.069;
mu = 0.011 * 1440; % convert from /min to /day to match lambdas
sat = Inf;
n0 = [10000000 0 0]; % 10^6 healthy cells
tspan = linspace(0, duration, 500);

[t,n] = ode45(@(t,n) simplified_model(t,n,a,b,@I_133_dose,total_dosage,mu,gamma,lambda, ...
    0,0,sat), tspan, n0);
end

function [t,n] = Dan_G_gen_uptake(lambda_s, lambda_d, total_dosage, duration)
a = 0.473;
b = 0.206;
lambda = 0.139;
gamma = 0.069;
mu = 0.011 * 1440; % convert from /min to /day to match lambdas
sat = Inf;
n0 = [10000000 0 0]; % 10^6 healthy cells
tspan = linspace(0, duration, 500);

[t,n] = ode45(@(t,n) model(t,n,a,b,@I_133_uptake,total_dosage,mu,gamma,lambda, ...
    lambda_s,lambda_d,sat), tspan, n0);
end

function [t,n] = Dan_G_simp_uptake(total_dosage, duration)
a = 0.473;
b = 0.206;
lambda = 0.139;
gamma = 0.069;
mu = 0.011 * 1440; % convert from /min to /day to match lambdas
sat = Inf;
n0 = [10000000 0 0]; % 10^6 healthy cells
tspan = linspace(0, duration, 500);

[t,n] = ode45(@(t,n) simplified_model(t,n,a,b,@I_133_uptake,total_dosage,mu,gamma,lambda, ...
    0,0,sat), tspan, n0);
end

%% Reduced model
function plot_reduced_expl(a, b)
% Plots the number of healthy cells and the number of sub-lethally damaged
% cells vs time using the reduced model for several values of a 
% fixed dose rate

% Set starting tumor size of 10^6 healthy cells
n0 = [10000000 0 0];

figure;
set(gcf,'units','normalized','outerposition',[0 0 0.7 0.3])

% r = 10
tspan = [0 1];
[t,n] = ode45(@(t,n) reduced_model(t,n,a*10,b*10), tspan, n0);
subplot(1, 4, 1); hold on
a1 = plot(t, n(:,1), "LineWidth", 2); M1 = "healthy cells";
a2 = plot(t, n(:,2), "LineWidth", 2); M2 = "sub-lethally damaged cells";
ylim([0 inf])
xlim([0 inf])
title("r = 10")
xlabel("time (days)")
ylabel("num cells")
legend([a1, a2], [M1, M2])

% r = 1
tspan = [0 12];
[t,n] = ode45(@(t,n) reduced_model(t,n,a,b), tspan, n0);
subplot(1, 4, 2); hold on
a1 = plot(t, n(:,1), "LineWidth", 2); M1 = "healthy cells";
a2 = plot(t, n(:,2), "LineWidth", 2); M2 = "sub-lethally damaged cells";
ylim([0 inf])
xlim([0 inf])
title("r = 1")
xlabel("time (days)")
ylabel("num cells")
legend([a1, a2], [M1, M2])


% r = 0.1
tspan = [0 90];
[t,n] = ode45(@(t,n) reduced_model(t,n,a*0.1,b*0.1), tspan, n0);
subplot(1, 4, 3); hold on
a1 = plot(t, n(:,1), "LineWidth", 2); M1 = "healthy cells";
a2 = plot(t, n(:,2), "LineWidth", 2); M2 = "sub-lethally damaged cells";
ylim([0 inf])
xlim([0 inf])
title("r = 0.1")
xlabel("time (days)")
ylabel("num cells")
legend([a1, a2], [M1, M2])

% r = 0.01
tspan = [0 900];
[t,n] = ode45(@(t,n) reduced_model(t,n,a*0.01,b*0.01), tspan, n0);
subplot(1, 4, 4); hold on
a1 = plot(t, n(:,1), "LineWidth", 2); M1 = "healthy cells";
a2 = plot(t, n(:,2), "LineWidth", 2); M2 = "sub-lethally damaged cells";
ylim([0 inf])
xlim([0 inf])
title("r = 0.01")
xlabel("time (days)")
ylabel("num cells")
legend([a1, a2], [M1, M2])

end


%% Dosage functions

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


%% Model functions


function nsd_dt = model(t, n, a, b, dosage, dose_amount, u, g, l, ls, ld, sat)
% The full model. Dosage is a function that computes the amount of
% radiation in the system at time t.

r = dosage(t, dose_amount);
nsd_dt = zeros(3, 1);
% [N, Ns, Nd]
sat_term = (1 - (n(1) + n(2) + n(3)) / sat);
nsd_dt(1) = l * n(1) * sat_term + u * n(2) - (a + b) * r * n(1);
nsd_dt(2) = ls * n(2) * sat_term + b * r * n(1) - u * n(2) - (a + b) * r * n(2);
nsd_dt(3) = ld * n(3) * sat_term + a * r * (n(1) + n(2)) + b * r * n(2) - g * n(3);
end


function nsd_dt = simplified_model(t, n, a, b, dosage, dose_amount, u, g, l, ls, ld, sat)
% The simplified model. 
% Dosage is a function that computes the amount of
% radiation in the system at time t.
% Sub-lethally damaged and doomed cells do not reproduce.

r = dosage(t, dose_amount);
nsd_dt = zeros(3, 1);
% [N, Ns, Nd]
sat_term = (1 - (n(1) + n(2) + n(3)) / sat);
nsd_dt(1) = l * n(1) * sat_term + u * n(2) - (a + b) * r * n(1);
nsd_dt(2) = b * r * n(1) - u * n(2) - (a + b) * r * n(2);
nsd_dt(3) = a * r * (n(1) + n(2)) + b * r * n(2) - g * n(3);
end


function nsd_dt = reduced_model(t,n,a,b)
% The reduced model.
% The effects of sub lethal repair damage are ignored, 
% as well as proliferation and elimination of doomed cells,
% and the dose rate is constant.

nsd_dt = zeros(3,1);
% [N, Ns, Nd]
nsd_dt(1) = -(a+b)*n(1);
nsd_dt(2) = b*n(1) - (a+b)*n(2);
nsd_dt(3) = a*n(1) + (a+b)*n(2);
end