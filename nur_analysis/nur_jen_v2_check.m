function R = nur_jen_v2_check(n_runs, param_set)
% NUR_JEN_V2_CHECK  Diagnostic of Jen's updated code (Drive: jen/jen_v2).
%
% Lives in nur_analysis/ inside a copy of jen_v2 and only CALLS Jen's
% functions in the parent folder (aDMDc, dmf_get_params, simulate_dmf,
% make_stimulus, ...). Nothing in her code is changed.
%
% What changed in jen_v2 compared with GitHub master:
%   - DMD on ALL 40 units (B_all, 40 x 304) with a new aDMDc.m function
%   - self-excitation only ON in areas where "lateral" is on (w_EE_on)
%   - old parameters (g_ff 1.0, g_fb 0.5, I0 0.3; master has 0.5/0.2/0.3255)
%   - main_exec: dt = 1, only scenario A
%
% This file runs three checks:
%
%   CHECK 1 - time step. Scenario A simulated with dt = 1 (as in jen_v2
%             main_exec) and dt = 0.01. With dt = 1 s the Euler step is 10x
%             the NMDA time constant (0.1 s), so S jumps and gets clipped at
%             0 or 1 instead of following the model.
%
%   CHECK 2 - the SVD step inside aDMDc. Reports how many singular values
%             of Omega are (numerically) zero, how many aDMDc keeps with its
%             current rule, and the smallest one it divides by.
%
%   CHECK 3 - recovery on random networks (independent tests, as in
%             nur_v3_random), using Jen's model rules. Every run gets random
%             feedforward / feedback strengths and random lateral on/off.
%             Three ways of estimating:
%               'Jen aDMDc (as is)'  : her aDMDc on all 40 units, then the
%                                      block means her main_exec reports
%               'units, SVD fixed'   : same 40-unit DMD, fixed rank rule
%               'areas, SVD fixed'   : V1 / V2 averages, fixed rank rule
%                                      (what my earlier diagnostics used)
%             Scores: Spearman rho true vs estimated for feedforward,
%             feedback, asymmetry (ff - fb) and leak (true fb vs estimated
%             ff, should be ~0).
%
% HOW TO RUN (in MATLAB):
%   cd ~/Documents/MATLAB/cn_pipeline_jen_v2/nur_analysis
%   nur_jen_v2_check                 % 40 random networks, Jen's parameters
%   nur_jen_v2_check(40, 'master')   % same, with master's I0 / w_EE
%
% Notes: simulations use dt = 0.01 (see check 1). BOLD uses the same
% Balloon-Windkessel equations as balloonWindkessel.m, integrated with
% Euler for speed (Jen's version uses ode45 with a waitbar per unit, which
% would open 40 waitbars per run). Her aDMDc builds a new random HGR basis
% every call, so rng(0) is set before each call to make it reproducible
% and identical to the basis used by the other two methods.

if nargin < 1 || isempty(n_runs),    n_runs = 40; end
if nargin < 2 || isempty(param_set), param_set = 'jen'; end

here   = fileparts(mfilename('fullpath'));
parent = fileparts(here);
fprintf('Running from: %s\n\n', here);
addpath(parent);
addpath(fullfile(parent, 'CNI_toolbox', 'code', 'matlab'));
assert(exist('aDMDc', 'file') == 2, 'aDMDc.m (jen_v2) not found in %s', parent);

p = dmf_get_params();
if strcmpi(param_set, 'master')          % values from GitHub master (Elena, 28 Sep)
    p.I0 = 0.3255;  p.w_EE_on = 0.3;
end
old = cd(parent);  setup_stimulus_grid;  cd(old);   % stimulus, A, X, Y, n_frames, T
n1 = 20;  n2 = 20;  area = [ones(1,n1), 2*ones(1,n2)];  N = numel(area);
p.prf = define_v1_retinotopy(ceil(sqrt(n1)), ceil(n1/ceil(sqrt(n1))), 10, 10, 1.5, 1);
p.prf = p.prf(1:n1);
fprintf('Parameters (%s): g_lat = %.2f, w_EE_on = %.2f, I0 = %.4f\n\n', ...
    param_set, p.g_lat, p.w_EE_on, p.I0);

%% ===================== CHECK 1: time step =============================
fprintf('=== CHECK 1: time step (scenario A) ===\n');
[C, w_EE] = make_scenario_connectivity('A', area, p);
dts = [1, 0.01];
figure('Name', 'Check 1: dt = 1 vs dt = 0.01', 'Position', [40 40 1100 450]);
for k = 1:2
    dt = dts(k);
    t_vec = (0:round(T/dt)-1) * dt;
    I_ext = make_stimulus(N, t_vec, p, A, X, Y);
    rng(1);
    [S_t, t] = simulate_dmf(C, w_EE, p, T, dt, [], I_ext);
    clipped = mean(S_t(:) <= 0 | S_t(:) >= 1) * 100;
    fprintf('dt = %-5g : mean S V1 = %.3f, V2 = %.3f | %5.1f%% of values clipped at 0 or 1\n', ...
        dt, mean(mean(S_t(area==1,:))), mean(mean(S_t(area==2,:))), clipped);
    subplot(1, 2, k);
    plot(t, mean(S_t(area==1,:),1), 'b', t, mean(S_t(area==2,:),1), 'r');
    xlabel('time (s)'); ylabel('S (area mean)'); ylim([0 1]);
    title(sprintf('dt = %g s  (%.0f%% clipped)', dt, clipped)); legend('V1','V2');
end
fprintf('-> use dt = 0.01 (or smaller); dt = 1 is 10x the NMDA time constant.\n\n');

%% ===================== CHECK 2: SVD step ===============================
fprintf('=== CHECK 2: SVD step inside aDMDc (scenario A, dt = 0.01) ===\n');
dt = 0.01;  t_vec = (0:round(T/dt)-1) * dt;
I_ext = make_stimulus(N, t_vec, p, A, X, Y);
rng(1);
[S_t, t] = simulate_dmf(C, w_EE, p, T, dt, [], I_ext);
B_all = bold_all(S_t, dt, p.frame_dt, n_frames);            % 40 x 304
gam = hgr_basis(p);
shash = (double(stimulus)' * gam)';
Om = [zscore(B_all(:,1:end-1)')'; zscore(shash(:,1:end-1)')'; zscore(shash(:,2:end)')'];
sv = svd(Om);
b_jen = 1;                                                   % size(Sig,2)/size(Sig,1) after econ SVD
thr_jen = optimal_SVHT_coef(b_jen, diag(sv));
kept = sv(sv >= thr_jen);
fprintf('Omega is %d x %d; %d of %d singular values are numerically zero (< 1e-8 x largest)\n', ...
    size(Om,1), size(Om,2), sum(sv < 1e-8*max(sv)), numel(sv));
fprintf('aDMDc as is keeps %d components; the smallest it divides by is %.2e\n', ...
    numel(kept), min(kept));
fprintf('  (median of all singular values = %.2e, so the threshold is ~0)\n', median(sv));
fprintf('fixed rule keeps %d components\n\n', fixed_rank(Om, sv));

%% ===================== CHECK 3: random networks ========================
fprintf('=== CHECK 3: recovery on %d random networks (TR = %g s) ===\n', n_runs, p.frame_dt);
methods = {'Jen aDMDc (as is)', 'units, SVD fixed', 'areas, SVD fixed'};
rows = {};
for run = 1:n_runs
    rng(5000 + run);
    g   = rand(1, 2);                     % ff, fb in [0, 1]
    g(rand(1, 2) < 0.2) = 0;              % sometimes absent
    lat = rand(1, 2) < 0.5;               % lateral on/off in V1, V2
    pr = p;  pr.g_ff = g(1);  pr.g_fb = g(2);
    [Cr, wr] = jen_connectivity([lat(1) 1; 1 lat(2)], area, pr);
    rng(1000 + run);
    S_t = simulate_dmf(Cr, wr, pr, T, dt, [], I_ext);
    B_all = bold_all(S_t, dt, p.frame_dt, n_frames);
    Xs = B_all(:, 1:end-1);  Xsp = B_all(:, 2:end);

    % (1) Jen's function, exactly as written
    rng(0);
    evalc('Aj = aDMDc(Xs, Xsp, stimulus, p);');          % evalc hides its disp(rtil)
    % (2) same 40-unit model, fixed rank
    Au = admdc_fixed(B_all, shash);
    % (3) area averages, fixed rank
    Ba = [mean(B_all(area==1,:),1); mean(B_all(area==2,:),1)];
    Aa = admdc_fixed(Ba, shash);

    est = [block(Aj, area, 2, 1), block(Aj, area, 1, 2); ...
           block(Au, area, 2, 1), block(Au, area, 1, 2); ...
           Aa(2,1), Aa(1,2)];
    for m = 1:3
        rows(end+1,:) = {methods{m}, run, g(1), est(m,1), g(2), est(m,2), lat(1), lat(2)}; %#ok<AGROW>
    end
    if mod(run, 10) == 0, fprintf('  %d / %d networks done\n', run, n_runs); end
end
R = cell2table(rows, 'VariableNames', {'method','run','g_ff','ff_est','g_fb','fb_est','lat_V1','lat_V2'});
writetable(R, fullfile(here, 'jen_v2_check_results.csv'));

RHO = zeros(3, 4);  PV = RHO;
rng(123);
for m = 1:3
    q = strcmp(R.method, methods{m});
    pairs = {R.g_ff(q), R.ff_est(q); R.g_fb(q), R.fb_est(q); ...
             R.g_ff(q) - R.g_fb(q), R.ff_est(q) - R.fb_est(q); R.g_fb(q), R.ff_est(q)};
    for k = 1:4
        [RHO(m,k), PV(m,k)] = spearman_perm(pairs{k,1}, pairs{k,2}, 2000);
    end
end
fprintf('\nSpearman rho (true vs estimated), * = p<0.05. leak should be ~0.\n');
fprintf('%-20s | %9s | %9s | %9s | %9s\n', 'method', 'ff', 'fb', 'asym', 'leak');
for m = 1:3
    fprintf('%-20s |', methods{m});
    for k = 1:4, fprintf(' %+5.2f %s  |', RHO(m,k), star(PV(m,k))); end
    fprintf('\n');
end

%% Figures
figure('Name', 'Check 3: recovery on random networks', 'Position', [60 60 1300 800]);
subplot(2, 3, 1);
bar(RHO(:, 1:3)');  set(gca, 'XTickLabel', {'ff', 'fb', 'asymmetry'});
ylim([-0.6 1.05]); ylabel('Spearman \rho'); legend(methods, 'Location', 'southoutside');
title('Recovery (higher = better)');
for m = 1:3
    q = strcmp(R.method, methods{m});
    subplot(2, 3, 3 + m);
    plot(R.g_ff(q), R.ff_est(q), 'bo', R.g_fb(q), R.fb_est(q), 'r^');
    xlabel('true strength'); ylabel('estimate');
    title(sprintf('%s\nff \\rho = %.2f, fb \\rho = %.2f', methods{m}, RHO(m,1), RHO(m,2)));
    legend('feedforward', 'feedback', 'Location', 'best');
end
end


%% =====================================================================
%  Local functions
%  =====================================================================

function [C, w_EE] = jen_connectivity(Mflag, area, p)
% Same code as Jen's make_scenario_connectivity.m, but with the on/off
% pattern (Mflag = [lat1 fb; ff lat2]) passed in directly, so random
% combinations can be used.
Gain = [p.g_lat p.g_fb; p.g_ff p.g_lat];
N = numel(area);
C = Mflag(area, area) .* Gain(area, area);
C(1:N+1:end) = 0;
for a = 1:2
    for b = 1:2
        ii = area == a;  jj = area == b;
        n  = nnz(jj) - (a == b);
        if n > 0, C(ii, jj) = C(ii, jj) / n; end
    end
end
lateral_on = diag(Mflag);
w_EE = p.w_EE_on * lateral_on(area);
end

function v = block(A, area, to, from)
% mean of the block of A from area 'from' into area 'to' (as in Jen's main_exec)
Ab = A(area == to, area == from);
v = mean(Ab(:));
end

function B = bold_all(S, dt, TR, n_frames)
% Balloon-Windkessel (same equations/parameters as balloonWindkessel.m) for
% every unit, then sampled every TR -> N x n_frames (like Jen's B_all)
kappa=0.65; gamma=0.41; tau_H=0.98; alpha=0.32; rho=0.34; V0=0.02;
n = size(S,1);  s = zeros(n,1); f = ones(n,1); v = ones(n,1); q = ones(n,1);
Bf = zeros(size(S));
for k = 1:size(S,2)
    Bf(:,k) = V0*(7*rho*(1-q) + 2*(1-q./v) + (2*rho-0.2)*(1-v));
    ds = S(:,k) - kappa*s - gamma*(f-1);
    dv = (f - v.^(1/alpha)) / tau_H;
    dq = ((f/rho).*(1-(1-rho).^(1./f)) - q.*v.^(1/alpha-1)) / tau_H;
    s = s + dt*ds;  f = max(f + dt*s, 1e-6);  v = v + dt*dv;  q = q + dt*dq;
end
B = Bf(:, 1:round(TR/dt):end);
B = B(:, 1:n_frames);
end

function gam = hgr_basis(p)
% same HGR settings as Jen's aDMDc, with rng(0)
rng(0);
hp.fwhm = 0.15; hp.eta = 0.1; hp.f_sampling = 1/p.frame_dt;
hp.r_stimulus = p.n_px; hp.n_features = 250; hp.n_gaussians = 5; hp.n_voxels = 2;
gam = HGR(hp).get_features;
end

function r = fixed_rank(Om, sv)
svk   = sv(sv > max(sv) * 1e-8);
beta  = min(size(Om)) / max(size(Om));
omega = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43;
r = sum(svk >= omega * median(svk));
end

function A = admdc_fixed(Bx, shash)
% Jen's aDMDc model (z-scored X, Xp, Ups, dUps), with the fixed rank rule
X  = zscore(Bx(:,1:end-1)')';  Xp = zscore(Bx(:,2:end)')';
Ctl = [shash(:,1:end-1); shash(:,2:end)];
Ctl = Ctl(std(Ctl, 0, 2) > 1e-12, :);
Om  = [X; zscore(Ctl')'];
[U, S, V] = svd(Om, 'econ');
sv = diag(S);
r  = min(max(fixed_rank(Om, sv), size(X,1)), sum(sv > max(sv)*1e-8));
A  = Xp * V(:,1:r) * diag(1 ./ sv(1:r)) * U(1:size(X,1), 1:r)';
end

function [rho, p] = spearman_perm(x, y, n_perm)
rx = ranks(x);  ry = ranks(y);
rho = pearson(rx, ry);
if isnan(rho), p = 1; return; end
null = zeros(n_perm, 1);
for b = 1:n_perm, null(b) = pearson(rx(randperm(numel(rx))), ry); end
p = (1 + sum(abs(null) >= abs(rho))) / (1 + n_perm);
end

function r = ranks(x)
x = x(:);  [~, i] = sort(x);  r = zeros(size(x));  r(i) = 1:numel(x);
u = unique(x);
for k = 1:numel(u)
    t = x == u(k);
    if sum(t) > 1, r(t) = mean(r(t)); end
end
end

function c = pearson(a, b)
a = a - mean(a);  b = b - mean(b);
den = sqrt(sum(a.^2) * sum(b.^2));
if den == 0, c = NaN; else, c = sum(a .* b) / den; end
end

function s = star(p)
if p < 0.05, s = '*'; else, s = ' '; end
end
