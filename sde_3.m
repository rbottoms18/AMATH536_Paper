% Ronan Bottoms
% AMATH 536
% June 5, 2024

delta_n = [1 0 0;
           0 1 0;
           0 0 1;
           1 -1 0;
           -1 1 0;
           -1 0 1;
           0 -1 1;
           0 0 -1];

l = 0.139;
a = 0.473;
b = 0.206;
u = 15.54;
g = 0.069;
sat = 1000000000000;

% Adjustable growth rate parameters
ls = 0.139;
ld = 5.55*10^(-3);
total_dosage = 10;
n = [1000000 0 0];

timesteps = 1000000;
delta_t = 1;
times = linspace(0, timesteps * delta_t, timesteps / 1000);
pop = zeros(timesteps / 1000, 3);

% Run simulation
for t=1:timesteps
    % Compute r
    r = I_133_uptake(total_dosage, t/1440);

    % Vector of probabilities
    p = [l * n(1) * (1 - (n(1) + n(2) + n(3))/sat) * delta_t;
         ls * n(2) * (1 - (n(1) + n(2) + n(3))/sat) * delta_t;
         ld * n(3) * (1 - (n(1) + n(2) + n(3))/sat) * delta_t;
         u * n(2) * delta_t;
         b * r * n(1) * delta_t;
         a * r * n(1) * delta_t;
         (a + b) * r * n(2) * delta_t;
         g * n(3) * delta_t];
    
    % Randomly draw which event occurs using a weighted random.
    s = sum(p);
    x = unifrnd(0, s);
    event_index = 8;
    for i = 1:8
        if (x < p(i))
            event_index = 1;
        end
        x = x - p(i);
    end

    n = n + delta_n(event_index);
    
    % Record every 1000th iteration
    if (mod(t, 1000) == 0)
        pop(t / 1000,:) = n;
    end
end

figure; hold on;
plot(times, pop(:, 1))
plot(times, pop(:, 2))
plot(times, pop(:, 3))


function dose_t = I_133_uptake(t, amount)
% Returns an array of the dosage at time t (days) for a single injection
% of amount 'amount' of iodine 133 that follows an uptake model.
dose_t = amount - amount * power(1/2, t / 8);
end


