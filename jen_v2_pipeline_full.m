%% JEN_V2_PIPELINE_FULL
% Jen's jen_v2 pipeline with all fixes and improvements, in one file.
%
%   stimulus --(pRFs)--> DMF (Gilson/Deco) --> Balloon-Windkessel BOLD
%            --> aDMDc (hashed stimulus, HGR features) --> connective
%            fields (A) and V1 pRFs (F), scored against the ground truth.
%
% Needs only stimulus.mat in the same folder (no CNI toolbox: the HGR
% feature generator is included below). Every function this script uses
% is a local function at the bottom of this file, so older versions of
% simulate_dmf.m etc. in the same folder are NOT used.
%
% MODES (run_mode, set in section 0):
%   'single' : one scenario, prints CF + pRF tables, plots
%   'all'    : scenarios A-F, ON vs OFF comparison per connection, saves
%              results_all_scenarios_full.mat, 3 CSVs and 2 figures
%   'sanity' : Option B check. Linear CF data (no DMF), like the Connective
%              Field Modelling chapter, to verify aDMDc itself; rank sweep
%
% CHANGES vs Jen's jen_v2 (see project notes jen_v2_pipeline_check.md)
%   1. dt = 0.1 ms (was 1 s -> Euler unstable with tau_NMDA = 0.1 s).
%      Output averaged into 10 ms bins; stimulus input stored per frame.
%   2. Hemodynamic lag: hashed stimulus convolved with the two-gamma HRF
%      before aDMDc (BOLD lags the stimulus by ~2 TRs).
%   3. SVHT rank cut-off: beta from the shape of Omega (was always 1).
%   4. aDMDc z-scores only inside the SVD; A, B, F returned in original
%      units (F back-projection now accounts for the feature scaling).
%   5. First 5 BOLD frames dropped (Balloon-Windkessel start-up transient).
%   6. Balloon-Windkessel: interpolant built once (was rebuilt every step).
%   7. Separate ranks: A at p.dmd_rank_cf (connective fields), F at
%      p.dmd_rank (pRFs). The SVHT rank (~136) recovers neither.
%   8. Octavian's topographic connectivity (GitHub 1c99b5d, 30 Sep):
%      Gaussian CFs over pRF distance; w_EE constant for all nodes.
%      GitHub parameter values (Elena 28 Sep, Octavian 30 Sep).
%   9. Scoring: per-node correlation / Jaccard similarity of A's rows vs.
%      the true weights (+ null), instead of block means; ON vs OFF across
%      scenarios. Reproducible seeds.

clear; close all; clc;

%% 0. Settings
run_mode = 'single';   % 'single' | 'all' | 'sanity'
scenario = 'A';        % used in 'single' mode
fast     = false;      % true -> dt = 1 ms (same results in tests, ~10x faster)
n1 = 20;  n2 = 20;     % V1 and V2 nodes

%% 1. Parameters, stimulus, visual-field grid
p = dmf_get_params();
if fast, p.dt = 1e-3; end
[stimulus, X, Y, n_frames, frame_t, T, p] = setup_stimulus_grid(p);

area = [ones(1, n1), 2*ones(1, n2)];
N    = numel(area);

%% 2. V1 pRFs, V2 anchors, stimulus input (shared by all scenarios)
n1_x = ceil(sqrt(n1));  n1_y = ceil(n1/n1_x);
prf  = define_v1_retinotopy(n1_x, n1_y, p.fov_radius, p.fov_radius, 1.5, 1);
rng(0);
keep   = sort(randperm(n1_x * n1_y, n1));   % n1 = 20 fills the 5x4 grid
p.prf  = prf(keep);
anchor = randperm(n1, n2);                  % V1 node each V2 node is centred on

% one column per frame (N x n_frames); simulate_dmf looks up the frame
I_ext = make_stimulus(N, frame_t, p, stimulus, X, Y);

% true V1 pRF images, for scoring F
true_prf = zeros(p.n_px^2, n1);
for v = 1:n1
    g = exp(-((X - p.prf(v).x0).^2 + (Y - p.prf(v).y0).^2) / (2*p.prf(v).sigma^2));
    true_prf(:,v) = g(:) / sum(g(:));
end

conn_names = {'V1->V2 (feedforward)','V2->V1 (feedback)', ...
              'V1->V1 (lateral)','V2->V2 (lateral)'};
here = fileparts(mfilename('fullpath'));

switch lower(run_mode)
%% ========================================================================
case 'single'
    res = run_one_scenario(scenario, p, area, n1, n2, N, anchor, I_ext, ...
                           T, n_frames, stimulus, true_prf);

    % block means (kept from Jen's version; A and C are in different
    % units, so only the sign/pattern is meaningful)
    C_full = res.C;  C_full(1:N+1:end) = p.w_EE;
    blocks = {'V1->V1 (lateral)', 1:n1, 1:n1;  'V2->V2 (lateral)', n1+1:N, n1+1:N;
              'V1->V2 (feedforward)', n1+1:N, 1:n1;  'V2->V1 (feedback)', 1:n1, n1+1:N};
    fprintf('\n%-22s %10s %10s\n', 'Block', 'mean(A)', 'mean(C)');
    for k = 1:size(blocks,1)
        rows = blocks{k,2}; cols = blocks{k,3};
        Ablk = res.A(rows,cols);  Cblk = C_full(rows,cols);
        if isequal(rows,cols), Ablk = Ablk(~eye(numel(rows))); Cblk = Cblk(~eye(numel(rows))); end
        fprintf('%-22s %10.4f %10.4f\n', blocks{k,1}, mean(Ablk(:)), mean(Cblk(:)));
    end

    print_cf_table(res.cf, conn_names, scenario, res.rank_cf);
    fprintf('V1 pRF recovery (rank %d): median r = %.3f (range %.3f to %.3f)\n', ...
        res.rank_prf, median(res.prf_r), min(res.prf_r), max(res.prf_r));

    figure('Position',[100 100 1100 450]);
    subplot(1,2,1); imagesc(res.W.ff); axis square; colorbar;
    title('True feedforward CFs (V2 x V1)'); xlabel('V1 node'); ylabel('V2 node');
    subplot(1,2,2); imagesc(res.A(n1+1:N,1:n1)); axis square; colorbar;
    title(sprintf('Estimated A, V1 \\rightarrow V2 (rank %d)', res.rank_cf));
    xlabel('V1 node'); ylabel('V2 node');

    figure('Position',[150 150 1100 450]);
    subplot(1,2,1); imagesc(res.A); axis square; colorbar; title('Estimated A');
    subplot(1,2,2); imagesc(C_full); axis square; colorbar; title('Ground truth (C + w_{EE})');

    figure;
    histogram(res.prf_r, 10);
    xlabel('correlation (true vs. estimated V1 pRF)'); ylabel('# nodes');
    title(sprintf('V1 pRF recovery, median r = %.2f', median(res.prf_r)));

%% ========================================================================
case 'all'
    scenarios = {'A','B','C','D','E','F'};
    nS = numel(scenarios);  nC = numel(conn_names);
    clear results
    for s = 1:nS
        fprintf('\n===== Scenario %s (%d/%d) =====\n', scenarios{s}, s, nS);
        results(s) = run_one_scenario(scenarios{s}, p, area, n1, n2, N, anchor, ...
                                      I_ext, T, n_frames, stimulus, true_prf); %#ok<SAGROW>
        print_cf_table(results(s).cf, conn_names, scenarios{s}, results(s).rank_cf);
        fprintf('V1 pRF median r = %.2f\n', median(results(s).prf_r));
    end

    % table: one row per scenario x connection
    Scenario = {}; Connection = {}; On = []; CF_r = []; Null_r = []; CF_JS = [];
    for s = 1:nS
        for k = 1:nC
            Scenario{end+1,1}   = scenarios{s};          %#ok<SAGROW>
            Connection{end+1,1} = conn_names{k};         %#ok<SAGROW>
            On(end+1,1)     = results(s).cf(k,1);        %#ok<SAGROW>
            CF_r(end+1,1)   = results(s).cf(k,2);        %#ok<SAGROW>
            Null_r(end+1,1) = results(s).cf(k,3);        %#ok<SAGROW>
            CF_JS(end+1,1)  = results(s).cf(k,4);        %#ok<SAGROW>
        end
    end
    T_scores = table(Scenario, Connection, On, CF_r, Null_r, CF_JS);

    % ON vs OFF per connection
    Mean_r_ON = zeros(nC,1); Mean_r_OFF = Mean_r_ON; Diff_r = Mean_r_ON;
    Mean_JS_ON = Mean_r_ON;  Mean_JS_OFF = Mean_r_ON;
    Scen_ON = cell(nC,1); Scen_OFF = cell(nC,1);
    for k = 1:nC
        rows = strcmp(T_scores.Connection, conn_names{k});
        on = rows & T_scores.On == 1;  off = rows & T_scores.On == 0;
        Mean_r_ON(k)  = mean(T_scores.CF_r(on));   Mean_r_OFF(k)  = mean(T_scores.CF_r(off));
        Mean_JS_ON(k) = mean(T_scores.CF_JS(on));  Mean_JS_OFF(k) = mean(T_scores.CF_JS(off));
        Diff_r(k) = Mean_r_ON(k) - Mean_r_OFF(k);
        Scen_ON{k}  = strjoin(T_scores.Scenario(on)', ',');
        Scen_OFF{k} = strjoin(T_scores.Scenario(off)', ',');
    end
    Connection = conn_names';
    T_summary = table(Connection, Scen_ON, Scen_OFF, Mean_r_ON, Mean_r_OFF, Diff_r, ...
                      Mean_JS_ON, Mean_JS_OFF);

    Scenario = scenarios';
    PRF_median_r = arrayfun(@(r) median(r.prf_r), results)';
    PRF_min_r    = arrayfun(@(r) min(r.prf_r),    results)';
    PRF_max_r    = arrayfun(@(r) max(r.prf_r),    results)';
    T_prf = table(Scenario, PRF_median_r, PRF_min_r, PRF_max_r);

    fprintf('\n===== ON vs OFF (CF rank %d) =====\n', p.dmd_rank_cf);  disp(T_summary);
    fprintf('===== V1 pRF recovery (rank %d) =====\n', p.dmd_rank);   disp(T_prf);

    f1 = figure('Name','A per scenario','Position',[50 50 1400 800]);
    for s = 1:nS
        subplot(2,3,s); imagesc(results(s).A); axis square; colorbar; hold on;
        xline(n1+0.5,'w-','LineWidth',1.5); yline(n1+0.5,'w-','LineWidth',1.5);
        fl = results(s).W.flags;
        title(sprintf('%s  ff%d fb%d latV1%d latV2%d', scenarios{s}, fl(2,1), fl(1,2), fl(1,1), fl(2,2)));
        xlabel('from node'); ylabel('to node');
    end
    f2 = figure('Name','CF score ON vs OFF','Position',[100 100 900 450]);
    bar(categorical(conn_names, conn_names), [Mean_r_ON Mean_r_OFF]);
    legend({'connection ON','connection OFF'}, 'Location','northoutside','Orientation','horizontal');
    ylabel(sprintf('mean CF r (rank %d)', p.dmd_rank_cf)); grid on;

    save(fullfile(here,'results_all_scenarios_full.mat'), 'results','T_scores', ...
         'T_summary','T_prf','p','anchor','scenarios','conn_names','-v7.3');
    writetable(T_scores,  fullfile(here,'scenario_scores_full.csv'));
    writetable(T_summary, fullfile(here,'scenario_onoff_summary_full.csv'));
    writetable(T_prf,     fullfile(here,'prf_scores_full.csv'));
    exportgraphics(f1, fullfile(here,'fig_A_blocks_full.png'), 'Resolution',150);
    exportgraphics(f2, fullfile(here,'fig_onoff_full.png'),    'Resolution',150);
    fprintf('\nSaved results_all_scenarios_full.mat, 3 CSVs and 2 PNGs in %s\n', here);

%% ========================================================================
case 'sanity'
    run_sanity_check(p, stimulus, X, Y, n_frames, true_prf, n1, n2);

otherwise
    error('run_mode must be ''single'', ''all'' or ''sanity''.');
end


%% ************************************************************************
%  LOCAL FUNCTIONS
%  ************************************************************************

%% ---------------- parameters & stimulus ---------------------------------
function p = dmf_get_params()
% DMF parameters (Gilson 2016, Table 3) + connectivity, stimulus, numerics
p.a        = 270;     % (n/C)
p.b        = 108;     % (Hz)
p.c        = 0.154;   % (s)   (called d in Deco et al. 2013)
p.beta     = 0.641;   % kinetic parameter
p.J        = 0.261;   % scaling constant
p.tau_NMDA = 100e-3;  % NMDA decay time constant (s)
p.I0       = 0.3255;  % background input (nA) [GitHub, Elena 28 Sep]
p.G        = 1.0;     % global coupling

% connectivity (GitHub: Elena 28 Sep, Octavian 30 Sep)
p.g_ff     = 0.7;     % feedforward strength  V1 -> V2
p.g_fb     = 0.4;     % feedback strength     V2 -> V1
p.g_lat    = 0.1;     % lateral (within-area) strength
p.w_EE     = 0.35;     % self-excitation, same for every node
p.s_ff     = 2;       % feedforward CF width (deg)
p.s_fb     = 3;       % feedback width (deg)
p.s_lat    = 1;       % lateral width (deg)

% noise (added to dS/dt, Euler-Maruyama)
p.sigma    = 0.01;

% external stimulus (Wong & Wang 2006; into V1 only)
p.J_ext    = 5.2e-4;  % AMPA coupling of external input (nA/Hz)
p.mu0      = 80;      % stimulus strength (Hz)

% numerics
p.dt       = 1e-4;    % 0.1 ms (Wong & Wang 2006); tau_NMDA = 100 ms
p.dt_store = 10e-3;   % outputs averaged into 10 ms bins
p.n_drop   = 5;       % BOLD frames dropped at the start (BW transient)

% aDMDc
p.stim_lag    = 'hrf';  % 'hrf' = two-gamma HRF, or integer = shift in TRs
p.dmd_rank    = 50;     % rank for pRFs (F); [] = SVHT
p.dmd_rank_cf = 10;     % rank for connective fields (A)
end

function [stimulus, X, Y, n_frames, frame_t, T, p] = setup_stimulus_grid(p)
% Bar stimulus (22500 x 304: 150x150 binary frames as columns) and the
% visual-field coordinates of every pixel (deg; x left->right, y top->bottom)
folder = fileparts(mfilename('fullpath'));
S = load(fullfile(folder, 'stimulus.mat'));
stimulus = double(S.stimulus);
p.n_px       = 150;
p.fov_radius = 10;     % image spans -10 to +10 deg
p.frame_dt   = 2;      % s per frame = TR
R = p.fov_radius;
[X, Y]   = meshgrid(linspace(-R, R, p.n_px), linspace(R, -R, p.n_px));
n_frames = size(stimulus, 2);
frame_t  = (0:n_frames-1) * p.frame_dt;
T        = n_frames * p.frame_dt;    % 608 s
end

function [prf, n_voxels] = define_v1_retinotopy(n_x, n_y, extent_x, extent_y, sigma, n_exp)
% V1 pRF centres on a Cartesian grid over the visual field
x_vals = linspace(-extent_x, extent_x, n_x);
y_vals = linspace(-extent_y, extent_y, n_y);
[Xg, Yg] = meshgrid(x_vals, y_vals);
x0 = Xg(:);  y0 = Yg(:);
n_voxels = numel(x0);
prf(n_voxels) = struct('x0', [], 'y0', [], 'sigma', [], 'n', []);
for v = 1:n_voxels
    prf(v).x0 = x0(v);  prf(v).y0 = y0(v);
    prf(v).sigma = sigma;  prf(v).n = n_exp;
end
end

function r = pRF_response(stim, X, Y, prf)
% Overlap between each frame and a normalised Gaussian pRF, in [0,1]
% (Dumoulin & Wandell 2008, eq. 2)
g = exp(-((X - prf.x0).^2 + (Y - prf.y0).^2) / (2*prf.sigma^2));
g = g(:) / sum(g(:));
r = g' * double(stim);                 % 1 x n_frames
end

function I_ext = make_stimulus(N, t_vec, p, stim, X, Y)
% External current into V1 (first numel(p.prf) nodes), J_ext*mu0*overlap
n_frames = size(stim, 2);
n_v1 = numel(p.prf);
idx = min(floor(t_vec / p.frame_dt) + 1, n_frames);
assert(idx(1) == 1 && idx(end) == n_frames, 'Simulation length does not match the stimulus length.');
I_ext = zeros(N, numel(t_vec));
for v = 1:n_v1
    r = pRF_response(stim, X, Y, p.prf(v));
    I_ext(v, :) = p.J_ext * p.mu0 * r(idx);
end
end

%% ---------------- connectivity (Octavian, GitHub 1c99b5d) ---------------
function [C, w_EE, W] = make_scenario_connectivity(scenario, p, n1, n2, N, anchor)
% Topographic ground-truth connectivity for scenarios A-F.
% C(i,j) = connection FROM node j INTO node i. V1 = 1:n1, V2 = n1+1:N.
%   Scenario   feedforward  feedback  lateral V1  lateral V2
%      A           x           x
%      B                                  x           x
%      C           x           x                      x
%      D           x                      x           x
%      E                       x          x           x
%      F           x           x          x
% W returns the unscaled weight matrices (ground-truth CFs) and flags.
prf_pos_v1 = [[p.prf.x0]', [p.prf.y0]'];
switch upper(scenario)   % Mflag(1,1) latV1, (1,2) fb, (2,1) ff, (2,2) latV2
    case 'A', Mflag = [0 1; 1 0];
    case 'B', Mflag = [1 0; 0 1];
    case 'C', Mflag = [0 1; 1 1];
    case 'D', Mflag = [1 0; 1 1];
    case 'E', Mflag = [1 1; 0 1];
    case 'F', Mflag = [1 1; 1 0];
    otherwise, error('Unknown scenario "%s". Use A-F.', scenario);
end
Gain = [p.g_lat p.g_fb; p.g_ff p.g_lat];
nrm  = @(M) M ./ sum(M, 2);
d2   = @(P, Q) (P(:,1) - Q(:,1)').^2 + (P(:,2) - Q(:,2)').^2;

% feedforward: each V2 node = Gaussian patch of V1 around its anchor
W_ff = nrm(exp(-d2(prf_pos_v1(anchor,:), prf_pos_v1) / (2*p.s_ff^2)));
prf_pos_v2 = W_ff * prf_pos_v1;                 % (larger) V2 pRF centres

% feedback: V1 node gets input from V2 nodes with nearby pRF centres
W_fb = nrm(exp(-d2(prf_pos_v1, prf_pos_v2) / (2*p.s_fb^2)));

% lateral: Gaussian over pRF distance within each area, no self-loops
W_lat_v1 = exp(-d2(prf_pos_v1, prf_pos_v1) / (2*p.s_lat^2));
W_lat_v1(1:n1+1:end) = 0;  W_lat_v1 = nrm(W_lat_v1);
W_lat_v2 = exp(-d2(prf_pos_v2, prf_pos_v2) / (2*p.s_lat^2));
W_lat_v2(1:n2+1:end) = 0;  W_lat_v2 = nrm(W_lat_v2);

i1 = 1:n1;  i2 = n1+1:N;
C = zeros(N);
C(i1,i1) = Mflag(1,1) * Gain(1,1) * W_lat_v1;
C(i1,i2) = Mflag(1,2) * Gain(1,2) * W_fb;
C(i2,i1) = Mflag(2,1) * Gain(2,1) * W_ff;
C(i2,i2) = Mflag(2,2) * Gain(2,2) * W_lat_v2;
w_EE = p.w_EE * ones(N,1);

W.ff = W_ff;  W.fb = W_fb;  W.lat_v1 = W_lat_v1;  W.lat_v2 = W_lat_v2;
W.flags = Mflag;  W.prf_pos_v2 = prf_pos_v2;
end

%% ---------------- DMF ----------------------------------------------------
function [S_t, t, u_t, r_t] = simulate_dmf(C, w_EE, p, T, dt, S0, I_ext, dt_store)
% Euler-Maruyama integration of the DMF (Deco 2013 / Gilson 2016).
% I_ext: (N x n_steps) or (N x n_frames, one column per frame_dt).
% Outputs averaged over dt_store windows.
N = size(C, 1);
n_steps = round(T / dt);
if isempty(S0), S0 = 0.1 * ones(N, 1); end
if nargin < 8 || isempty(dt_store), dt_store = dt; end
every   = round(dt_store / dt);
n_store = floor(n_steps / every);
t = (0:n_store-1) * every * dt;

per_step = size(I_ext, 2) == n_steps;
n_cols   = size(I_ext, 2);

S_t = zeros(N, n_store);  u_t = S_t;  r_t = S_t;
S = S0;  w_EE = w_EE(:);
accS = zeros(N,1); accU = accS; accR = accS;
noise_sd = p.sigma * sqrt(dt);
for k = 1:n_steps
    if per_step
        Ik = I_ext(:, k);
    else
        Ik = I_ext(:, min(floor((k-1)*dt / p.frame_dt) + 1, n_cols));
    end
    [dSdt, u, h] = dmf_activity_change(S, C, w_EE, p, Ik);
    S = S + dSdt * dt + noise_sd * randn(N, 1);
    S = min(max(S, 0), 1);
    accS = accS + S;  accU = accU + u;  accR = accR + h;
    if mod(k, every) == 0
        j = k / every;
        S_t(:, j) = accS / every;  u_t(:, j) = accU / every;  r_t(:, j) = accR / every;
        accS(:) = 0;  accU(:) = 0;  accR(:) = 0;
    end
end
end

function [dS_dt, u, h] = dmf_activity_change(S, C, w_EE, p, I_ext)
% dS/dt (eq. 22a)
u = dmf_input(S, C, w_EE, p, I_ext);
h = dmf_activation(u, p);
dS_dt = -S ./ p.tau_NMDA + p.beta .* (1 - S) .* h;
end

function u = dmf_input(S, C, w_EE, p, I_ext)
% total input current (eq. 22b)
if nargin < 5 || isempty(I_ext), I_ext = zeros(size(S)); end
u = p.J .* (w_EE .* S + p.G .* (C * S)) + p.I0 + I_ext;
end

function h = dmf_activation(x, p)
% sigmoidal transfer function H(x) (eq. 22)
num = p.a .* x - p.b;
den = 1 - exp(-p.c .* num);
h = num ./ den;
h(den == 0) = 1 / p.c;     % 0/0 limit
end

%% ---------------- Balloon-Windkessel -------------------------------------
function B = balloonWindkessel(S, t)
% BOLD from synaptic activity (eq. 23, Table 3 parameters), ode45
kappa = 0.65; gamma_ = 0.41; tau_H = 0.98; alpha = 0.32; rho = 0.34; V0 = 0.02;
y0 = [0; 1; 1; 1];                                  % s, f, v, q at rest
S_of_t = griddedInterpolant(t(:), S(:), 'linear', 'nearest');   % built once
model = @(tc, y) [S_of_t(tc) - kappa*y(1) - gamma_*(y(2) - 1);
                  y(1);
                  (y(2) - y(3)^(1/alpha)) / tau_H;
                  (1/tau_H) * ((y(2)/rho) * (1 - (1-rho)^(1/y(2))) - y(4) * y(3)^(1/alpha - 1))];
[~, ySol] = ode45(model, t, y0);
v = ySol(:,3);  q = ySol(:,4);
B = V0 .* (7*rho .* (1 - q) + 2 .* (1 - q ./ v) + (2*rho - 0.2) .* (1 - v));
end

%% ---------------- aDMDc --------------------------------------------------
function [A, B, F, info] = aDMDc(X, Xp, stimulus, p)
% Algebraic DMD with control:  X(:,k+1) = A X(:,k) + B u(:,k) + F u(:,k+1)
% u = stimulus hashed into HGR features, delayed to match the BOLD lag.
% A (n x n), B (n x q), F (n x n_pixels) in original units.
gam   = hgr_features(p.n_px, 250, 5, 0.15);       % n_pixels x 250
shash = double(stimulus)' * gam;                  % n_frames x 250

% hemodynamic lag
if ischar(p.stim_lag) || isstring(p.stim_lag)
    t_hrf = (0:p.frame_dt:34-1)';
    shash = filter(two_gamma(t_hrf), 1, shash);   % causal convolution in time
else
    L = round(p.stim_lag);
    shash = [zeros(L, size(shash,2)); shash(1:end-L,:)];
end
shash = shash';                                   % 250 x n_frames
n_snap = size(X,2) + 1;
shash  = shash(:, end-n_snap+1:end);              % same frames as the BOLD
U  = shash(:,1:end-1);
dU = shash(:,2:end);

% z-score rows for the SVD only
[Xz,  ~, sX ] = zscore(X,  0, 2);
[Xpz, ~, sXp] = zscore(Xp, 0, 2);
[Uz,  ~, sU ] = zscore(U,  0, 2);
[dUz, ~, sdU] = zscore(dU, 0, 2);
sU(sU == 0) = 1;  sdU(sdU == 0) = 1;

Omega = [Xz; Uz; dUz];
[Uo, Sig, V] = svd(Omega, 'econ');
sv = diag(Sig);

% rank cut-off (SVHT, beta from Omega's shape) or forced rank
beta = min(size(Omega)) / max(size(Omega));
threshold = optimal_SVHT_coef(beta) * median(sv);
svht_rank = sum(sv >= threshold);
if isfield(p,'dmd_rank') && ~isempty(p.dmd_rank), rtil = p.dmd_rank; else, rtil = svht_rank; end
fprintf('aDMDc: Omega %dx%d, beta = %.3f, SVHT rank = %d, used rank = %d\n', ...
    size(Omega,1), size(Omega,2), beta, svht_rank, rtil);

Util = Uo(:,1:rtil);  Sigtil = Sig(1:rtil,1:rtil);  Vtil = V(:,1:rtil);
n = size(X,1);  q = size(U,1);
K  = Xpz * Vtil / Sigtil;
Az = K * Util(1:n,:)';
Bz = K * Util(n+1:n+q,:)';
Fz = K * Util(n+q+1:n+2*q,:)';

% back to original units
A = sXp .* Az ./ sX';
B = sXp .* Bz ./ sU';
F = sXp .* Fz ./ sdU';
F = (gam * F')';                                  % n x n_pixels

info.rank = rtil;  info.svht_rank = svht_rank;
info.threshold = threshold;  info.beta = beta;  info.sv = sv;
end

function gam = hgr_features(r_stimulus, n_features, n_gaussians, fwhm)
% Hashed Gaussian features, identical to HGR.create_gamma (CNI toolbox,
% Bhat et al. 2021): each feature = sum of n_gaussians Gaussians at random
% pixel positions, normalised to sum 1. Returns n_pixels x n_features.
Yg = linspace(0, r_stimulus, r_stimulus)' * ones(1, r_stimulus);
Xg = ones(r_stimulus, 1) * linspace(0, r_stimulus, r_stimulus);
sigma = fwhm * r_stimulus / (2 * sqrt(2 * log(2)));
pix_id = linspace(0, r_stimulus^2 - 1, n_features * n_gaussians);
pix_id = pix_id(randperm(n_features * n_gaussians));
x = floor(pix_id / r_stimulus) + 1;
y = mod(pix_id, r_stimulus) + 1;
gam = zeros(r_stimulus^2, n_features);
for i = 0:n_features-1
    g = zeros(r_stimulus);
    for j = 1:n_gaussians
        k = i*n_gaussians + j;
        g = g + exp(-((Xg - x(k)).^2 + (Yg - y(k)).^2) ./ (2 * sigma^2));
    end
    gam(:, i+1) = g(:) / sum(g(:));
end
end

function h = two_gamma(t)
% canonical two-gamma HRF (same as HGR.two_gamma)
h = (6*t.^5.*exp(-t))./gamma(6) - 1/6*(16*t.^15.*exp(-t))/gamma(16);
end

function omega = optimal_SVHT_coef(beta)
% Gavish & Donoho (2014) optimal hard threshold coefficient, noise level
% unknown: threshold = omega(beta) * median(singular values)
w = (8 * beta) ./ (beta + 1 + sqrt(beta.^2 + 14 * beta + 1));
lambda_star = sqrt(2 * (beta + 1) + w);
omega = lambda_star ./ sqrt(median_marcenko_pastur(beta));
end

function med = median_marcenko_pastur(beta)
lobnd = (1 - sqrt(beta))^2;  hibnd = (1 + sqrt(beta))^2;
topSpec = hibnd;  botSpec = lobnd;
dens = @(x) sqrt(max((topSpec - x).*(x - botSpec), 0)) ./ (beta .* x) ./ (2*pi);
MarPas = @(x0) 1 - integral(dens, x0, topSpec);
change = 1;
while change && (hibnd - lobnd > .001)
    change = 0;
    xs = linspace(lobnd, hibnd, 5);
    ys = arrayfun(MarPas, xs);
    if any(ys < 0.5), lobnd = max(xs(ys < 0.5)); change = 1; end
    if any(ys > 0.5), hibnd = min(xs(ys > 0.5)); change = 1; end
end
med = (hibnd + lobnd) / 2;
end

%% ---------------- one full scenario --------------------------------------
function res = run_one_scenario(sc, p, area, n1, n2, N, anchor, I_ext, T, n_frames, stimulus, true_prf) %#ok<INUSL>
% DMF -> BOLD -> aDMDc -> scores, for one scenario. Same noise seed and
% HGR features in every scenario.
[C, w_EE, W] = make_scenario_connectivity(sc, p, n1, n2, N, anchor);

rng(1);
tic;
[S_t, t] = simulate_dmf(C, w_EE, p, T, p.dt, [], I_ext, p.dt_store);
fprintf('Scenario %s: DMF done in %.0f s\n', sc, toc);

ds = round(p.frame_dt / p.dt_store);
B_all = zeros(N, n_frames);
for v = 1:N
    Bv = balloonWindkessel(S_t(v,:), t);
    B_all(v,:) = Bv(1:ds:end);
    if mod(v, 10) == 0, fprintf('  BOLD %d/%d nodes\n', v, N); end
end
B_all = B_all(:, p.n_drop+1:end);
Xs = B_all(:,1:end-1);  Xsp = B_all(:,2:end);

p_cf = p;  p_cf.dmd_rank = p.dmd_rank_cf;
rng(100); [A_cf, ~, ~, info_cf] = aDMDc(Xs, Xsp, stimulus, p_cf);   % CFs
rng(100); [~, ~, F, info_prf]   = aDMDc(Xs, Xsp, stimulus, p);      % pRFs

prf_r = zeros(n1,1);
for v = 1:n1
    prf_r(v) = corr(F(v,:)', true_prf(:,v));
end

res.scenario = sc;  res.C = C;  res.w_EE = w_EE;  res.W = W;
res.BOLD = B_all;   res.A = A_cf;  res.F = F;
res.rank_cf = info_cf.rank;  res.rank_prf = info_prf.rank;
res.cf = score_connections(A_cf, W, n1, N);
res.prf_r = prf_r;
end

%% ---------------- scoring ------------------------------------------------
function tab = score_connections(A, W, n1, N)
% one row per connection [ff; fb; latV1; latV2]:
% [on, median CF r, median null r, median Jaccard]. Null = same rows
% matched to the wrong receiving nodes. Compares SHAPES (A and C are in
% different units).
blk = {n1+1:N, 1:n1,   W.ff,     W.flags(2,1);
       1:n1,   n1+1:N, W.fb,     W.flags(1,2);
       1:n1,   1:n1,   W.lat_v1, W.flags(1,1);
       n1+1:N, n1+1:N, W.lat_v2, W.flags(2,2)};
tab = zeros(size(blk,1), 4);
for k = 1:size(blk,1)
    rows = blk{k,1};  cols = blk{k,2};  Wt = blk{k,3};
    Ablk = A(rows, cols);
    if isequal(rows, cols), Ablk(1:numel(rows)+1:end) = 0; end
    Wnull = circshift(Wt, 7, 1);
    r = zeros(numel(rows),1);  nl = r;  js = r;
    for i = 1:numel(rows)
        r(i)  = corr(Ablk(i,:)', Wt(i,:)');
        nl(i) = corr(Ablk(i,:)', Wnull(i,:)');
        js(i) = jaccard_sim(Ablk(i,:), Wt(i,:));
    end
    tab(k,:) = [blk{k,4}, median(r), median(nl), median(js)];
end
end

function js = jaccard_sim(est, tru)
% chapter's Jaccard similarity: sum(min)/sum(max) after clipping negative
% weights to 0 and scaling both vectors to a peak of 1
est = max(est, 0);  tru = max(tru, 0);
if max(est) == 0 || max(tru) == 0, js = 0; return; end
est = est / max(est);  tru = tru / max(tru);
js = sum(min(est, tru)) / sum(max(est, tru));
end

function print_cf_table(tab, conn_names, sc, rank_cf)
fprintf('\nScenario %s, CF rank %d\n', sc, rank_cf);
fprintf('%-22s %4s %8s %8s %8s\n', 'Connection', 'on?', 'CF r', 'null r', 'CF JS');
for k = 1:numel(conn_names)
    fprintf('%-22s %4d %8.2f %8.2f %8.2f\n', conn_names{k}, tab(k,:));
end
end

%% ---------------- Option B sanity check ----------------------------------
function run_sanity_check(p, stimulus, X, Y, n_frames, true_prf, n1, n2)
% Data built like the Connective Field Modelling chapter (sec. 3.2.3):
% V1 = pRF overlap, V2 = W * V1 in the same TR (Gaussian CF), BOLD = HRF
% * neural + noise. If aDMDc recovers W and the pRFs here, the method is
% fine. Rank sweep shows the CF vs pRF trade-off.
rng(1);
N = n1 + n2;
sigma_cf  = 5;       % CF width (deg)
noise_rel = 0.1;     % noise std as fraction of BOLD std
ranks     = [NaN 10 20 30 50 80];   % NaN = SVHT

r_v1 = zeros(n1, n_frames);
for v = 1:n1, r_v1(v,:) = pRF_response(stimulus, X, Y, p.prf(v)); end
x0 = [p.prf.x0]';  y0 = [p.prf.y0]';
W = exp(-((x0 - x0').^2 + (y0 - y0').^2) / (2*sigma_cf^2));
W = W ./ sum(W, 2);
neural = [r_v1; W * r_v1];

BOLD = filter(two_gamma((0:p.frame_dt:34-1)'), 1, neural, [], 2);
BOLD = BOLD + noise_rel * std(BOLD(:)) * randn(size(BOLD));
BOLD = BOLD(:, p.n_drop+1:end);
Xs = BOLD(:,1:end-1);  Xsp = BOLD(:,2:end);

W_null = circshift(W, 7, 1);
fprintf('\n%-10s %8s %10s %8s %8s\n', 'rank', 'CF r', 'CF null r', 'CF JS', 'pRF r');
best = -Inf;
for k = 1:numel(ranks)
    if isnan(ranks(k)), p.dmd_rank = []; else, p.dmd_rank = ranks(k); end
    rng(100);
    [Ahat, ~, F, info] = aDMDc(Xs, Xsp, stimulus, p);
    A_cf = Ahat(n1+1:N, 1:n1);
    cf_r = zeros(n2,1); cf_n = cf_r; cf_js = cf_r;
    for j = 1:n2
        cf_r(j)  = corr(A_cf(j,:)', W(j,:)');
        cf_n(j)  = corr(A_cf(j,:)', W_null(j,:)');
        cf_js(j) = jaccard_sim(A_cf(j,:), W(j,:));
    end
    prf_r = zeros(n1,1);
    for v = 1:n1, prf_r(v) = corr(F(v,:)', true_prf(:,v)); end
    lbl = sprintf('%d', info.rank);  if isnan(ranks(k)), lbl = [lbl ' (SVHT)']; end
    fprintf('%-10s %8.2f %10.2f %8.2f %8.2f\n', lbl, median(cf_r), median(cf_n), ...
        median(cf_js), median(prf_r));
    if median(cf_r) > best, best = median(cf_r); A_best = A_cf; r_best = info.rank; end
end

figure('Name', sprintf('Sanity check: CFs, rank %d', r_best));
subplot(1,2,1); imagesc(W); axis square; colorbar; title('True CF weights W (V2 x V1)');
xlabel('V1 node'); ylabel('V2 node');
subplot(1,2,2); imagesc(A_best); axis square; colorbar; title('Estimated A (V1 \rightarrow V2)');
xlabel('V1 node'); ylabel('V2 node');
end
