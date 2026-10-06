function run_sanity_check(p, stimulus, X, Y, n_frames, true_prf, n1, n2, anchor)
% RUN_SANITY_CHECK  Option B check of aDMDc on linear CF data (no DMF).
% Data built like the Connective Field Modelling chapter (sec. 3.2.3):
% V1 = pRF overlap, V2 = W * V1 in the same TR (Gaussian CF around each
% V2 node's anchor in V1), BOLD = HRF * neural + noise. If aDMDc recovers
% W and the pRFs here, the method is fine. Rank sweep shows the CF vs pRF
% trade-off.
%
%   stimulus : (n_pixels x n_frames) the frames that were simulated
%   anchor   : (1 x n2) V1 node each V2 node is centred on

rng(1);
N = n1 + n2;
sigma_cf  = 5;       % CF width (deg)
noise_rel = 0.1;     % noise std as fraction of BOLD std
ranks     = [NaN 10 20 30 50 80];   % NaN = SVHT

%% Neural signals
r_v1 = zeros(n1, n_frames);
for v = 1:n1
    r_v1(v,:) = pRF_response(stimulus, X, Y, p.prf(v));
end

% CF weights (n2 x n1): each V2 node = Gaussian patch of V1 around its
% anchor, same construction as the feedforward connections in
% make_scenario_connectivity
prf_pos_v1 = [[p.prf.x0]', [p.prf.y0]'];
anchor_pos = prf_pos_v1(anchor, :);
dx = anchor_pos(:,1) - prf_pos_v1(:,1)';
dy = anchor_pos(:,2) - prf_pos_v1(:,2)';
W  = exp(-(dx.^2 + dy.^2) / (2*sigma_cf^2));
W  = W ./ sum(W, 2);

neural = [r_v1; W * r_v1];                        % N x n_frames

%% BOLD
BOLD = filter(two_gamma((0:p.frame_dt:34-1)'), 1, neural, [], 2);
BOLD = BOLD + noise_rel * std(BOLD(:)) * randn(size(BOLD));
BOLD = BOLD(:, p.n_drop+1:end);
Xs  = BOLD(:, 1:end-1);
Xsp = BOLD(:, 2:end);

%% Rank sweep
W_null = circshift(W, 7, 1);
fprintf('\n%-10s %8s %10s %8s %8s\n', 'rank', 'CF r', 'CF null r', 'CF JS', 'pRF r');
best = -Inf;
for k = 1:numel(ranks)
    if isnan(ranks(k)), p.dmd_rank = []; else, p.dmd_rank = ranks(k); end
    rng(100);
    [Ahat, ~, F, info] = aDMDc(Xs, Xsp, stimulus, p);
    A_cf = Ahat(n1+1:N, 1:n1);

    cf_r = zeros(n2, 1);
    cf_n = cf_r;
    cf_js = cf_r;
    for j = 1:n2
        cf_r(j)  = corr(A_cf(j,:)', W(j,:)');
        cf_n(j)  = corr(A_cf(j,:)', W_null(j,:)');
        cf_js(j) = jaccard_sim(A_cf(j,:), W(j,:));
    end

    prf_r = zeros(n1, 1);
    for v = 1:n1
        prf_r(v) = corr(F(v,:)', true_prf(:,v));
    end

    lbl = sprintf('%d', info.rank);
    if isnan(ranks(k)), lbl = [lbl ' (SVHT)']; end
    fprintf('%-10s %8.2f %10.2f %8.2f %8.2f\n', lbl, median(cf_r), median(cf_n), ...
        median(cf_js), median(prf_r));

    if median(cf_r) > best
        best   = median(cf_r);
        A_best = A_cf;
        r_best = info.rank;
    end
end

%% Plot best CF recovery
figure('Name', sprintf('Sanity check: CFs, rank %d', r_best));
subplot(1,2,1); imagesc(W); axis square; colorbar; title('True CF weights W (V2 x V1)');
xlabel('V1 node'); ylabel('V2 node');
subplot(1,2,2); imagesc(A_best); axis square; colorbar; title('Estimated A (V1 \rightarrow V2)');
xlabel('V1 node'); ylabel('V2 node');
end
