%% example_run_gilson_dmf.m
%
% Example driver for the Gilson/Deco DMF model
%
%   stimulus --(pRFs)--> DMF --> Balloon-Windkessel BOLD --> aDMDc
%            --> connective fields (A) and V1 pRFs (F), scored against
%                the ground truth
%
% Requires: dmf_get_params.m, setup_stimulus_grid.m, stimulus.mat,
%           define_V1_retinotopy.m, make_stimulus.m,
%           make_scenario_connectivity.m, pRF_response.m 
%           dmf_activation.m, dmf_input.m,
%           dmf_activity_change.m, simulate_dmf.m, balloonWindkessel.m,
%           run_one_scenario.m, aDMDc.m, hgr_features.m, two_gamma.m,
%           optimal_SVHT_coef.m, score_connections.m, jaccard_sim.m,
%           summarize_scenarios.m, run_sanity_check.m
%           (all in the same folder / path)
%
% MODES (run_mode, section 0):
%   'single' : one scenario, CF + pRF tables and plots
%   'all'    : scenarios A-F, ON vs OFF comparison per connection; saves
%              a .mat and 3 CSVs
%   All outputs go to output/ (data/, tables/, figures/); every open
%   figure is saved as output/figures/<run_mode>_<figure name>.png
%   'sanity' : aDMDc on linear CF data (no DMF), rank sweep

clear; close all; clc;

%% 0. Run mode
run_mode = 'all';   % 'single' | 'all' | 'sanity'
scenario = 'A';     % used in 'single' mode

%% 1. Load DMF parameters
p = dmf_get_params();

%% 2. Stimulus and visual-field grid
setup_stimulus_grid;       % creates A, X, Y, n_frames, frame_t, T

%% 3. Simulation settings
n1 = 100;  n2 = 50;
area = [ones(1, n1), 2*ones(1, n2)];
N    = numel(area);
iV1 = find(area == 1, 1); % this is always 1
iV2 = find(area == 2, 1); % this is always n1 + 1
dt       = 1e-3;
dt_store = 10e-3;   % simulate_dmf averages its outputs over bins of this length (s)
n_steps  = round(T / dt);
t_vec    = (0:n_steps-1) * dt;
S0       = [];


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
figure('Name', 'V1 pRFs'); imagesc(X(1,:), Y(:,1), G); axis xy equal tight; colorbar; hold on
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

% True V1 pRF images (same normalised Gaussian as in pRF_response), one
% column per V1 node -- ground truth for the aDMDc pRF estimates F
true_prf = zeros(p.n_px^2, n1);
for v = 1:n1
    g = exp(-((X - p.prf(v).x0).^2 + (Y - p.prf(v).y0).^2) / (2*p.prf(v).sigma^2));
    true_prf(:, v) = g(:) / sum(g(:));
end

% One column per frame (N x n_frames); simulate_dmf looks up the frame
I_ext = make_stimulus(N, frame_t, p, A, X, Y);

conn_names = {'V1->V2 (feedforward)', 'V2->V1 (feedback)', ...
              'V1->V1 (lateral)', 'V2->V2 (lateral)'};
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end

% output folders (created if missing)
out_dir   = fullfile(here, 'output');
data_dir  = fullfile(out_dir, 'data');      % .mat files
table_dir = fullfile(out_dir, 'tables');    % CSVs
fig_dir   = fullfile(out_dir, 'figures');   % PNGs
for d = {out_dir, data_dir, table_dir, fig_dir}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end

%% 5. Define structural connectivity (ground truth C) and run scenarios

% Running each scenario from the phd paper in the order found in figure 3.1
% Each scenario and the connectivity mtx is set up in make_scenario_connectivity
rng(0);
anchor = randperm(n1,n2);

if strcmpi(run_mode, 'sanity')
    run_sanity_check(p, A, X, Y, n_frames, true_prf, n1, n2, anchor);
    save_all_figures(fig_dir, run_mode);
    return
end

switch lower(run_mode)
    case 'all',    scenarios = {'A', 'B', 'C', 'D', 'E', 'F'};
    case 'single', scenarios = {scenario};
    otherwise,     error('run_mode must be ''single'', ''all'' or ''sanity''.');
end
nS = numel(scenarios);

clear results
for k = 1:nS
    sc = scenarios{k};
    fprintf('\n===== Scenario %s (%d/%d) =====\n', sc, k, nS);

    % seed k: same noise per scenario across runs (reproducible)
    results(k) = run_one_scenario(sc, p, n1, n2, N, anchor, I_ext, T, dt, ...
                     dt_store, S0, n_frames, A, true_prf, k); %#ok<SAGROW>

    fprintf('Scenario %s: mean S  V1 = %.3f, V2 = %.3f\n', ...
            sc, mean(results(k).S_t(iV1, :)), mean(results(k).S_t(iV2, :)));
    print_cf_table(results(k).cf, conn_names, sc, results(k).rank_cf);
    fprintf('V1 pRF recovery (rank %d): median r = %.2f (range %.2f to %.2f)\n', ...
            results(k).rank_prf, median(results(k).prf_r), ...
            min(results(k).prf_r), max(results(k).prf_r));
end


%% 5. Basic sanity checks / summary stats
fprintf('\nSimulation complete.\n');
fprintf('  N regions      : %d\n', N);
fprintf('  Duration       : %.1f s\n', T);
fprintf('  Time step      : %.4f s\n', dt);
fprintf('  Mean S (final) : %.4f\n', mean(results(end).S_t(:, end)));
fprintf('  Std  S (final) : %.4f\n', std(results(end).S_t(:, end)));

%% 6. Plot neural activity trajectories


figure('Name', 'S trajectories', 'Position', [100 100 1000 700]);
for k = 1:nS
    subplot(3, 2, k); hold on;
    plot(results(k).t, results(k).S_t(iV1, :), 'b');
    plot(results(k).t, results(k).S_t(iV2, :), 'r');
    xlim([0 T]);
    title(sprintf('Scenario %s', results(k).scenario));
    xlabel('Time (s)'); ylabel('S');
    if k == 1, legend({'V1', 'V2'}, 'Location', 'northeast'); end
end

%% 7. Plot firing rates H(u)
% This is just for one node in V1 and one node in V2
figure('Name', 'Firing rates', 'Position', [100 100 1000 700]);
for k = 1:nS
    subplot(3, 2, k); hold on;

    % shade the stimulus-on blocks
    yl = [0, max(results(k).r(:)) * 1.1];

    plot(results(k).t, results(k).r(iV1, :), 'b');
    plot(results(k).t, results(k).r(iV2, :), 'r');

    xlim([0 T]); ylim(yl);
    title(sprintf('Scenario %s', results(k).scenario));
    xlabel('Time (s)'); ylabel('Firing rate H(u) (Hz)');
    if k == 1, legend({'V1', 'V2'}, 'Location', 'northeast'); end
end

%% 8. Plot firing rates for every node (heatmap)
figure('Name', 'Firing rate heatmap', 'Position', [100 100 1000 700]);
for k = 1:nS
    subplot(3, 2, k);
    imagesc(results(k).t, 1:N, results(k).r);
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

%% 9. aDMDc results
switch lower(run_mode)
case 'single'
    res = results(1);

    % ground truth with w_EE on the diagonal
    C_full = res.C;
    C_full(1:N+1:end) = res.w_EE;

    % block means (A and C are in different units, so only the
    % sign/pattern is meaningful)
    blocks = {'V1->V1 (lateral)',     1:n1,   1:n1;
              'V2->V2 (lateral)',     n1+1:N, n1+1:N;
              'V1->V2 (feedforward)', n1+1:N, 1:n1;
              'V2->V1 (feedback)',    1:n1,   n1+1:N};
    fprintf('\n%-22s %10s %10s\n', 'Block', 'mean(A)', 'mean(C)');
    for k = 1:size(blocks,1)
        rows = blocks{k,2};
        cols = blocks{k,3};
        Ablk = res.A_est(rows, cols);
        Cblk = C_full(rows, cols);
        if isequal(rows, cols)
            Ablk = Ablk(~eye(numel(rows)));
            Cblk = Cblk(~eye(numel(rows)));
        end
        fprintf('%-22s %10.4f %10.4f\n', blocks{k,1}, mean(Ablk(:)), mean(Cblk(:)));
    end

    figure('Name', 'Feedforward CFs true vs estimated', 'Position', [100 100 1100 450]);
    subplot(1,2,1); imagesc(res.W.ff); axis square; colorbar;
    title('True feedforward CFs (V2 x V1)'); xlabel('V1 node'); ylabel('V2 node');
    subplot(1,2,2); imagesc(res.A_est(n1+1:N, 1:n1)); axis square; colorbar;
    title(sprintf('Estimated A, V1 \\rightarrow V2 (rank %d)', res.rank_cf));
    xlabel('V1 node'); ylabel('V2 node');

    figure('Name', 'A vs ground truth', 'Position', [150 150 1100 450]);
    subplot(1,2,1); imagesc(res.A_est); axis square; colorbar; title('Estimated A');
    subplot(1,2,2); imagesc(C_full);    axis square; colorbar; title('Ground truth (C + w_{EE})');

    figure('Name', 'pRF recovery');
    histogram(res.prf_r, 10);
    xlabel('correlation (true vs. estimated V1 pRF)'); ylabel('# nodes');
    title(sprintf('V1 pRF recovery, median r = %.2f', median(res.prf_r)));

case 'all'
    [T_scores, T_summary, T_prf] = summarize_scenarios(results, scenarios, conn_names);
    fprintf('\n===== ON vs OFF (CF rank %d) =====\n', p.dmd_rank_cf);  disp(T_summary);
    fprintf('===== V1 pRF recovery (rank %d) =====\n', p.dmd_rank);   disp(T_prf);

    % estimated A per scenario
    f1 = figure('Name', 'A per scenario', 'Position', [50 50 1400 800]);
    for s = 1:nS
        subplot(2,3,s); imagesc(results(s).A_est); axis square; colorbar; hold on;
        xline(n1+0.5, 'w-', 'LineWidth', 1.5); yline(n1+0.5, 'w-', 'LineWidth', 1.5);
        fl = results(s).W.flags;
        title(sprintf('%s  ff%d fb%d latV1%d latV2%d', scenarios{s}, fl(2,1), fl(1,2), fl(1,1), fl(2,2)));
        xlabel('from node'); ylabel('to node');
    end

    % ground truth (C + w_EE on the diagonal), same layout as f1, one
    % colour scale for all six panels so ON/OFF blocks are comparable
    C_gt = cell(1, nS);
    for s = 1:nS
        C_gt{s} = results(s).C;
        C_gt{s}(1:N+1:end) = results(s).w_EE;
    end
    gt_lim = [min(cellfun(@(M) min(M(:)), C_gt)), max(cellfun(@(M) max(M(:)), C_gt))];
    f3 = figure('Name', 'Ground truth per scenario', 'Position', [75 75 1400 800]);
    for s = 1:nS
        subplot(2,3,s); imagesc(C_gt{s}); axis square; colorbar; hold on;
        set(gca, 'CLim', gt_lim);
        xline(n1+0.5, 'w-', 'LineWidth', 1.5); yline(n1+0.5, 'w-', 'LineWidth', 1.5);
        fl = results(s).W.flags;
        title(sprintf('%s (truth)  ff%d fb%d latV1%d latV2%d', scenarios{s}, fl(2,1), fl(1,2), fl(1,1), fl(2,2)));
        xlabel('from node'); ylabel('to node');
    end

    % CF score, connection ON vs OFF
    f2 = figure('Name', 'CF score ON vs OFF', 'Position', [100 100 900 450]);
    bar(categorical(conn_names, conn_names), [T_summary.Mean_r_ON T_summary.Mean_r_OFF]);
    legend({'connection ON', 'connection OFF'}, 'Location', 'northoutside', 'Orientation', 'horizontal');
    ylabel(sprintf('mean CF r (rank %d)', p.dmd_rank_cf)); grid on;

    % save (S_t, u, r left out of the .mat: hundreds of MB for 150 nodes)
    out.results    = rmfield(results, {'S_t', 'u', 'r'});
    out.T_scores   = T_scores;
    out.T_summary  = T_summary;
    out.T_prf      = T_prf;
    out.p          = p;
    out.anchor     = anchor;
    out.scenarios  = scenarios;
    out.conn_names = conn_names;
    out.dt         = dt;
    out.dt_store   = dt_store;
    out.T          = T;
    save(fullfile(data_dir, 'dmf_six_scenarios.mat'), '-struct', 'out', '-v7.3');
    writetable(T_scores,  fullfile(table_dir, 'scenario_scores.csv'));
    writetable(T_summary, fullfile(table_dir, 'scenario_onoff_summary.csv'));
    writetable(T_prf,     fullfile(table_dir, 'prf_scores.csv'));
    fprintf('\nSaved dmf_six_scenarios.mat in %s and 3 CSVs in %s\n', data_dir, table_dir);
end

%% 10. Save every open figure
save_all_figures(fig_dir, run_mode);


%% Local functions
function print_cf_table(tab, conn_names, sc, rank_cf)
% one line per connection: on?, median CF r, median null r, median Jaccard
fprintf('\nScenario %s, CF rank %d\n', sc, rank_cf);
fprintf('%-22s %4s %8s %8s %8s\n', 'Connection', 'on?', 'CF r', 'null r', 'CF JS');
for k = 1:numel(conn_names)
    fprintf('%-22s %4d %8.2f %8.2f %8.2f\n', conn_names{k}, tab(k,:));
end
end

function save_all_figures(fig_dir, tag)
% saves every open figure as <tag>_<figure name>.png in fig_dir
figs = findall(groot, 'Type', 'figure');
for f = figs'
    name = f.Name;
    if isempty(name), name = sprintf('figure%d', f.Number); end
    name = regexprep(name, '[^\w-]+', '_');      % safe file name
    exportgraphics(f, fullfile(fig_dir, sprintf('%s_%s.png', tag, name)), 'Resolution', 150);
end
fprintf('Saved %d figures in %s\n', numel(figs), fig_dir);
end
