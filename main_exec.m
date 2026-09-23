%% example_run_gilson_dmf.m
%
% Example driver for the Gilson/Deco DMF model
%
% Requires: gilson_dmf_params.m, dmf_activation.m, dmf_input.m,
%           dmf_rhs.m, simulate_dmf.m  (all in the same folder / path)
% Test push git
clear; close all; clc;

%% 1. Load DMF parameters
p = dmf_get_params();

%% 2. Stimulus and visual-field grid
setup_stimulus_grid;       % creates A, X, Y, n_frames, T; adds fields to p

%% 3. Simulation settings
n1 = 20;  n2 = 20;
area = [ones(1, n1), 2*ones(1, n2)];
N    = numel(area);
iV1 = find(area == 1, 1);
iV2 = find(area == 2, 1);
dt      = 1e-3;
n_steps = round(T / dt);
t_vec   = (0:n_steps-1) * dt;
S0      = [];


%% 4. V1 pRF and stimulus input
prf = define_v1_retinotopy(5, 4, 10, 10, 1.5, 1);   % 40 V1 voxels, uniform sigma = 1.5 deg
p.prf = prf;

I_ext  = make_stimulus(N, t_vec, p, A, X, Y);

%% 5. Define structural connectivity (ground truth C)
% Replace this with curated connectivity matrix.
% Here: a small random sparse matrix, just for demonstration.

% Running each scenario from the phd paper in the order found in figure 3.1
% Each scenario and the connectivity mtx is set up in make_scenario_connectivity
scenarios = {'A', 'B', 'C', 'D', 'E', 'F'};
results = struct();

for k = 1:numel(scenarios)
    sc = scenarios{k};
    [C, w_EE] = make_scenario_connectivity(sc, area, p);

    rng(k);   % same noise per scenario across runs (reproducible)
    [S_t, t, u_t,r_t] = simulate_dmf(C, w_EE, p, T, dt, S0, I_ext);

    results(k).scenario = sc;
    results(k).C        = C;
    results(k).w_EE     = w_EE;
    results(k).S_t      = S_t;
    results(k).t        = t;
    results(k).I_ext = I_ext;
    results(k).u = u_t;
    results(k).r = r_t;

    fprintf('Scenario %s: mean S  V1 = %.3f, V2 = %.3f\n', ...
            sc, mean(S_t(iV1, :)), mean(S_t(iV2, :)));
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
    plot(t, results(k).S_t(iV1, :), 'b');
    plot(t, results(k).S_t(iV2, :), 'r');
    xlim([0 T]);
    title(sprintf('Scenario %s', results(k).scenario));
    xlabel('Time (s)'); ylabel('S');
    if k == 1, legend({'V1', 'V2'}, 'Location', 'northeast'); end
end

save('dmf_six_scenarios.mat', 'results', 'p', 'dt', 'T');

%% 7. Plot firing rates H(u)
figure('Position', [100 100 1000 700]);
for k = 1:numel(scenarios)
    subplot(3, 2, k); hold on;

    % shade the stimulus-on blocks
    yl = [0, max(results(k).r(:)) * 1.1];

    plot(t, results(k).r(iV1, :), 'b');
    plot(t, results(k).r(iV2, :), 'r');

    xlim([0 T]); ylim(yl);
    title(sprintf('Scenario %s', results(k).scenario));
    xlabel('Time (s)'); ylabel('Firing rate H(u) (Hz)');
    if k == 1, legend({'V1', 'V2'}, 'Location', 'northeast'); end
end