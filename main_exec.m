%% example_run_gilson_dmf.m
%
% Example driver for the Gilson/Deco DMF model
%
% Requires: dmf_get_params.m, setup_stimulus_grid.m, stimulus.mat,
%           define_V1_retinotopy.m, make_stimulus.m,
%           make_scenario_connectivity.m, pRF_response.m 
%           dmf_activation.m, dmf_input.m,
%           dmf_activity_change.m,
%           simulate_dmf.m  (all in the same folder / path)

clear; close all; clc;

%% 1. Load DMF parameters
p = dmf_get_params();

%% 2. Stimulus and visual-field grid
setup_stimulus_grid;       % creates A, X, Y, n_frames, T; adds fields to p

%% 3. Simulation settings
n1 = 100;  n2 = 50;
area = [ones(1, n1), 2*ones(1, n2)];
N    = numel(area);
iV1 = find(area == 1, 1); % this is always 1
iV2 = find(area == 2, 1); % this is always n1 + 1
dt      = 1e-3;
n_steps = round(T / dt);
t_vec   = (0:n_steps-1) * dt;
S0      = [];


%% 4. V1 pRF and stimulus input
% using the smallest grid that fits n1 points
n1_x = ceil(sqrt(n1));
n1_y = ceil(n1/n1_x);



prf = define_v1_retinotopy(n1_x, n1_y, p.fov_radius,p.fov_radius, 1.5, 1);   % n1 voxels, uniform sigma = 1.5 deg


% here to handle cases where n1 doesn't factor cleanly, randomizes gaps in
% the grid rather than having all of them in one spot (at the end)
%keep = sort(randperm(n1_x * n1_y, n1));
p.prf = prf;


% Visualize kept V1 pRFs, labeled by node index
G = zeros(size(X));
for v = 1:n1
    G = G + exp(-((X-p.prf(v).x0).^2 + (Y-p.prf(v).y0).^2)/(2*p.prf(v).sigma^2));
end
figure; imagesc(X(1,:), Y(:,1), G); axis xy equal tight; colorbar; hold on
th = linspace(0, 2*pi, 100);
for v = 1:n1
    x0 = p.prf(v).x0;  y0 = p.prf(v).y0;  s = p.prf(v).sigma;
    plot(x0 + s*cos(th), y0 + s*sin(th), 'w-');          % 1-sigma outline
    plot(x0, y0, 'w.', 'MarkerSize', 8);                 % centre
    text(x0, y0 + 0.5, sprintf('%d', v), 'Color', 'k', ...
        'BackgroundColor', [1 1 1 0.7], 'Margin', 1, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
xlabel('x (deg)'); ylabel('y (deg)'); title('V1 pRFs (label = node index)');

I_ext  = make_stimulus(N, t_vec, p, A, X, Y);

%% 5. Define structural connectivity (ground truth C)

% Running each scenario from the phd paper in the order found in figure 3.1
% Each scenario and the connectivity mtx is set up in make_scenario_connectivity
scenarios = {'A', 'B', 'C', 'D', 'E', 'F'};
results = struct();

rng(0);
anchor = randperm(n1,n2);

for k = 1:numel(scenarios)
    sc = scenarios{k};
    [C, w_EE] = make_scenario_connectivity(sc, area, p, n1, n2, N, anchor);

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

%save('dmf_six_scenarios.mat', 'results', 'p', 'dt', 'T');

%% 7. Plot firing rates H(u)
% This is just for one node in V1 and one node in V2
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

%% 8. Plot firing rates for every node (heatmap)
figure('Position', [100 100 1000 700]);
for k = 1:numel(scenarios)
    subplot(3, 2, k);
    imagesc(t, 1:N, results(k).r);
    set(gca, 'YDir', 'normal');
    colormap(gca, 'hot');
    cb = colorbar;
    cb.Label.String = 'Firing rate H(u) (Hz)';

    hold on;
    yline(n1 + 0.5, 'c-', 'LineWidth', 1.5);   % V1/V2 boundary

    xlim([0 T]);
    title(sprintf('Scenario %s', results(k).scenario));
    xlabel('Time (s)'); ylabel('Node index');
end