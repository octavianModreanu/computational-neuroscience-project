%% visualize_and_compare.m
%
% Run this AFTER main_exec.m (needs `results`, `p`, `T`, `dt`, `S0`,
% `stimulus` in the workspace) and after setup_stimulus_grid.m has
% loaded `X`, `Y`.
%
% Builds the HGR hashing basis once, runs the DMF->BOLD->aDMDc pipeline
% for every scenario, plots a per-scenario dashboard, and finishes with
% a summary of how well A recovers C across all six scenarios.
 
addpath("CNI_toolbox/code/matlab/")

%% Build the hashing basis ONCE (expensive, identical across scenarios)
hgr_params.fwhm        = 0.15;
hgr_params.eta         = 0.1;
hgr_params.f_sampling  = 1/p.frame_dt;
hgr_params.r_stimulus  = p.n_px;
hgr_params.n_features  = 250;
hgr_params.n_gaussians = 5;
hgr_params.n_voxels    = 2;
prf_mapper = HGR(hgr_params);
gam = prf_mapper.get_features;

%% Ground-truth V1 pRF image, for comparison against the estimated one
gt_prf_img = reshape( ...
    exp( -((X(:)-p.prf(1).x0).^2 + (Y(:)-p.prf(1).y0).^2) / (2*p.prf(1).sigma^2) ), ...
    p.n_px, p.n_px);

%% Run every scenario
scenarios = {results.scenario};
all_out = cell(1, numel(scenarios));
all_metrics = cell(1, numel(scenarios));

for k = 1:numel(scenarios)
    out = run_scenario_pipeline(k, results, p, T, dt, S0, stimulus, gam);
    all_out{k} = out;
    all_metrics{k} = compare_connectivity(out.C, out.A, out.scenario);

    %% Per-scenario dashboard
    figure('Name', ['Scenario ' out.scenario], 'Position', [50 50 1100 700]);

    subplot(2,3,1);
    plot(out.t, out.S_t(1,:), 'b', out.t, out.S_t(2,:), 'r');
    title('Synaptic activity S'); xlabel('t (s)'); legend('V1','V2');

    subplot(2,3,2);
    plot(out.t, out.r_t(1,:), 'b', out.t, out.r_t(2,:), 'r');
    title('Firing rate H(u)'); xlabel('t (s)'); legend('V1','V2');

    subplot(2,3,3);
    tr_axis = (0:numel(out.B1)-1) * p.frame_dt;
    plot(tr_axis, out.B1, 'b', tr_axis, out.B2, 'r');
    title('BOLD (downsampled to TR)'); xlabel('t (s)'); legend('V1','V2');

    subplot(2,3,4);
    imagesc(out.C); axis square; colorbar; caxis([-1 1]);
    title('Ground truth C'); set(gca,'XTick',1:2,'YTick',1:2, ...
        'XTickLabel',{'V1','V2'},'YTickLabel',{'V1','V2'});

    subplot(2,3,5);
    imagesc(out.A); axis square; colorbar; caxis([-1 1]);
    title(sprintf('Estimated A (fit R^2=%.2f)', out.fit_r2));
    set(gca,'XTick',1:2,'YTick',1:2,'XTickLabel',{'V1','V2'},'YTickLabel',{'V1','V2'});

    subplot(2,3,6);
    est_prf_img = reshape(out.F(1,:), p.n_px, p.n_px);
    imagesc(est_prf_img); axis square off; colormap(gca, 'hot');
    title('Estimated V1 pRF (pixel space)');

    sgtitle(sprintf('Scenario %s', out.scenario));
end

%% Cross-scenario summary: how well does A recover C?
ff_true = cellfun(@(m) m.ff_true, all_metrics);
ff_est  = cellfun(@(m) m.ff_est,  all_metrics);
fb_true = cellfun(@(m) m.fb_true, all_metrics);
fb_est  = cellfun(@(m) m.fb_est,  all_metrics);
fit_r2  = cellfun(@(o) o.fit_r2,  all_out);

figure('Name', 'Cross-scenario connectivity recovery', 'Position', [100 100 900 400]);

subplot(1,2,1);
bar([ff_true; ff_est]');
set(gca, 'XTickLabel', scenarios);
legend('true (V1->V2)', 'estimated', 'Location', 'best');
title('Feedforward strength'); ylabel('connectivity weight');

subplot(1,2,2);
bar([fb_true; fb_est]');
set(gca, 'XTickLabel', scenarios);
legend('true (V2->V1)', 'estimated', 'Location', 'best');
title('Feedback strength'); ylabel('connectivity weight');

sgtitle('True vs. estimated connectivity across scenarios A-F');

fprintf('\nModel fit (fraction of Xp variance explained) per scenario:\n');
for k = 1:numel(scenarios)
    fprintf('  %s: R^2 = %.3f\n', scenarios{k}, fit_r2(k));
end
