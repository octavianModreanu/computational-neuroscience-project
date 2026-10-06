function res = run_one_scenario(sc, p, n1, n2, N, anchor, I_ext, T, dt, dt_store, S0, n_frames, stimulus, true_prf, seed)
% RUN_ONE_SCENARIO  DMF -> Balloon-Windkessel BOLD -> aDMDc -> scores,
% for one scenario.
%
%   sc       : scenario letter 'A' ... 'F'
%   I_ext    : (N x n_frames) external input, one column per frame
%   T, dt    : simulation length and integration step (s)
%   dt_store : simulate_dmf output bin length (s)
%   S0       : initial condition ([] = default)
%   stimulus : (n_pixels x n_frames) the frames that were simulated
%   true_prf : (n_pixels x n1) true V1 pRF images, for scoring F
%   seed     : rng seed for the DMF noise
%
%   Returns res with the ground truth (C, w_EE, W), the simulation
%   (S_t, t, u, r, BOLD), the aDMDc estimates (A_est, F) and the scores.

[C, w_EE, W] = make_scenario_connectivity(sc, p, n1, n2, N, anchor);

%% DMF
rng(seed);
tic;
[S_t, t, u_t, r_t] = simulate_dmf(C, w_EE, p, T, dt, S0, I_ext, dt_store);
fprintf('Scenario %s: DMF done in %.0f s\n', sc, toc);

%% BOLD, sampled once per frame (TR)
ds = round(p.frame_dt / dt_store);
B_all = zeros(N, n_frames);
for v = 1:N
    Bv = balloonWindkessel(S_t(v,:), t);
    B_all(v,:) = Bv(1:ds:end);
    if mod(v, 25) == 0, fprintf('  BOLD %d/%d nodes\n', v, N); end
end
B_all = B_all(:, p.n_drop+1:end);      % drop Balloon-Windkessel start-up
Xs  = B_all(:, 1:end-1);
Xsp = B_all(:, 2:end);

%% aDMDc: A at the CF rank, F at the pRF rank (same HGR features: rng(100))
p_cf = p;
p_cf.dmd_rank = p.dmd_rank_cf;
rng(100); [A_cf, ~, ~, info_cf] = aDMDc(Xs, Xsp, stimulus, p_cf);   % CFs
rng(100); [~, ~, F, info_prf]   = aDMDc(Xs, Xsp, stimulus, p);      % pRFs

%% pRF recovery (V1)
prf_r = zeros(n1, 1);
for v = 1:n1
    prf_r(v) = corr(F(v,:)', true_prf(:,v));
end

%% Output
res.scenario = sc;
res.C        = C;
res.w_EE     = w_EE;
res.W        = W;
res.I_ext    = I_ext;
res.S_t      = S_t;
res.t        = t;
res.u        = u_t;
res.r        = r_t;
res.BOLD     = B_all;
res.A_est    = A_cf;
res.F        = F;
res.rank_cf  = info_cf.rank;
res.rank_prf = info_prf.rank;
res.cf       = score_connections(A_cf, W, n1, N);
res.prf_r    = prf_r;
end
