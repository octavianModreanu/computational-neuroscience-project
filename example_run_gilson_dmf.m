%% example_run_gilson_dmf.m
%
% Example driver for the Gilson/Deco DMF model
%
% Requires: gilson_dmf_params.m, dmf_activation.m, dmf_input.m,
%           dmf_rhs.m, simulate_dmf.m  (all in the same folder / path)

clear; close all; clc;

%% 1. Load DMF parameters
p = dmf_get_params();

%% 2. Simulation settings
N  = 2;
T  = 120;       % total time (s)
dt = 1e-3;      % integration step (s)
n_steps = round(T / dt);
t_vec = (0:n_steps-1) * dt; %time vector for stimulus (s)

S0 = [];      % use default initial condition (0.1 for all regions)

%% 3. Define structural connectivity (ground truth C)
% Replace this with curated connectivity matrix.
% Here: a small random sparse matrix, just for demonstration.

% Running each scenario from the phd paper in the order found in figure 3.1
% Each scenario and the connectivity mtx is set up in make_scenario_connectivity
scenarios = {'A', 'B', 'C', 'D', 'E', 'F'};
results = struct();

for k = 1:numel(scenarios)
    sc = scenarios{k};
    [C, w_EE] = make_scenario_connectivity(sc, p);

    rng(k);   % same noise per scenario across runs (reproducible)
    I_ext = make_stimulus(C,t_vec, p);
    [S_t, t] = simulate_dmf(C, w_EE, p, T, dt, [], I_ext);

    results(k).scenario = sc;
    results(k).C        = C;
    results(k).w_EE     = w_EE;
    results(k).S_t      = S_t;
    results(k).t        = t;
    results(k).I_ext = I_ext;

    fprintf('Scenario %s: mean S  V1 = %.3f, V2 = %.3f\n', ...
            sc, mean(S_t(1, :)), mean(S_t(2, :)));
end


%% 5. Basic sanity checks / summary stats
fprintf('Simulation complete.\n');
fprintf('  N regions      : %d\n', N);
fprintf('  Duration       : %.1f s\n', T);
fprintf('  Time step      : %.4f s\n', dt);
fprintf('  Mean S (final) : %.4f\n', mean(S_t(:, end)));
fprintf('  Std  S (final) : %.4f\n', std(S_t(:, end)));

%% 6. Plot neural activity trajectories
figure('Position', [100 100 1000 700]);
for k = 1:numel(scenarios)
    subplot(3, 2, k); hold on;
    plot(t, results(k).S_t(1, :), 'b');
    plot(t, results(k).S_t(2, :), 'r');
    xlim([0 T]);
    title(sprintf('Scenario %s', results(k).scenario));
    xlabel('Time (s)'); ylabel('S');
    if k == 1, legend({'V1', 'V2'}, 'Location', 'northeast'); end
end

save('dmf_six_scenarios.mat', 'results', 'p', 'dt', 'T');