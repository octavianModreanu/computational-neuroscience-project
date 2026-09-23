function out = run_scenario_pipeline(k, results, p, T, dt, S0, stimulus, gam)
% RUN_SCENARIO_PIPELINE  DMF -> BOLD -> hashed-stimulus aDMDc for one scenario.
%
%   k        : index into results (1..6, matches scenarios{k})
%   results  : struct array from main_exec.m (needs .C, .w_EE, .I_ext, .scenario)
%   p        : parameter struct (needs p.frame_dt, p.n_px)
%   T, dt    : simulation duration / integration step used in simulate_dmf
%   S0       : initial condition for simulate_dmf ([] = default)
%   stimulus : 22500 x n_frames raw bar movie (pixel space)
%   gam      : n_pixels x n_features hashing basis from HGR (build ONCE
%              outside this function/loop -- it's expensive and identical
%              across scenarios, since it only depends on the stimulus grid)
%
%   Returns out with fields:
%     .scenario, .C, .S_t, .t, .r_t, .B1, .B2, .A, .F, .Xpred, .fit_r2

    sc = results(k).scenario;
    fprintf('--- Scenario %s ---\n', sc);

    %% DMF activity
    [S_t, t, ~, r_t] = simulate_dmf(results(k).C, results(k).w_EE, p, T, dt, S0, results(k).I_ext);

    %% BOLD (label the progress bar so you can tell V1 from V2 while it runs)
    B1 = balloonWindkessel(S_t(1,:), t, sprintf('%s: V1 BOLD', sc));
    B2 = balloonWindkessel(S_t(2,:), t, sprintf('%s: V2 BOLD', sc));

    %% Downsample to TR so BOLD lines up with stimulus frames
    ds = round(p.frame_dt / dt);
    B1 = B1(1:ds:end);
    B2 = B2(1:ds:end);

    n_frames = min(numel(B1), size(stimulus, 2));
    B1 = B1(1:n_frames);
    B2 = B2(1:n_frames);

    data   = [B1(:), B2(:)];
    data_t = data';

    X  = data_t(:, 1:end-1);
    Xp = data_t(:, 2:end);

    %% Hashed stimulus as control input (pixel space -> feature space)
    shash = (stimulus(:, 1:n_frames)' * gam)';
    Ups   = zscore(shash(:, 1:end-1)')';
    dUps  = zscore(shash(:, 2:end)')';

    Omega = [X; Ups; dUps];

    %% aDMDc
    [U, Sig, V] = svd(Omega, 'econ');
    beta = size(Sig,2) / size(Sig,1);
    threshold = optimal_SVHT_coef(beta, Sig);
    rtil = sum(diag(Sig) >= threshold);

    Util = U(:,1:rtil); Sigtil = Sig(1:rtil,1:rtil); Vtil = V(:,1:rtil);

    n = size(X,1); q = size(Ups,1);
    U_1 = Util(1:n,:);
    U_2 = Util(n+1:n+q,:);
    U_3 = Util(n+q+1:n+q+q,:);

    A = Xp * Vtil * inv(Sigtil) * U_1';
    B = Xp * Vtil * inv(Sigtil) * U_2';
    F = Xp * Vtil * inv(Sigtil) * U_3';

    Xpred = A*X + B*Ups + F*dUps;

    % simple goodness-of-fit: fraction of variance in Xp explained by Xpred
    resid = Xp - Xpred;
    fit_r2 = 1 - sum(resid(:).^2) / sum((Xp(:) - mean(Xp(:))).^2);

    F = (gam * F')';   % back to pixel space, 2 x n_pixels

    out.scenario = sc;
    out.C        = results(k).C;
    out.S_t = S_t; out.t = t; out.r_t = r_t;
    out.B1 = B1; out.B2 = B2;
    out.A = A; out.F = F; out.Xpred = Xpred;
    out.fit_r2 = fit_r2;
end
