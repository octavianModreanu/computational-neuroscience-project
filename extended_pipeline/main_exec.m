%% i am trying my best 

clear; close all; clc;

%% 1. Load DMF parameters
p = dmf_get_params();

%% 2. Stimulus and visual-field grid
setup_stimulus_grid;       % creates A, X, Y, n_frames, T; adds fields to p

%% 3. Simulation settings
n1 = 20;  n2 = 20;
area = [ones(1, n1), 2*ones(1, n2)];
N    = numel(area);
iV1 = find(area == 1, 1); % this is always 1
iV2 = find(area == 2, 1); % this is always n1 + 1
dt      = 1;
n_steps = round(T / dt);
t_vec   = (0:n_steps-1) * dt;
S0      = [];

%% 4. V1 pRF and stimulus input
% using the smallest grid that fits n1 points
n1_x = ceil(sqrt(n1));
n1_y = ceil(n1/n1_x);

prf = define_v1_retinotopy(n1_x, n1_y, 10, 10, 1.5, 1);   % n1 voxels, uniform sigma = 1.5 deg
% here to handle cases where n1 doesn't factor cleanly, randomizes gaps in
% the grid rather than having all of them in one spot (at the end)
keep = sort(randperm(n1_x * n1_y, n1));
p.prf = prf(keep);

I_ext  = make_stimulus(N, t_vec, p, A, X, Y);

%% 5. Define structural connectivity (ground truth C)
% I chose to only use one scenario here, but my hope is that I can easily
% replace this with the original loop. 
 
[C, w_EE] = make_scenario_connectivity('A',area,p);

%% 6. Input structural connectivity into synaptic activity eq 22
[S_t, t, u_t,r_t] = simulate_dmf(C, w_EE, p, T, dt, S0, I_ext);

%% 7. Balloon-Windkessel for every node, then DOWNSAMPLE to TR
ds = round(p.frame_dt / dt);
B_all = zeros(N, n_frames);   % 40 x 304

for v = 1:N
    Bv = balloonWindkessel(S_t(v,:), t, sprintf('BW node %d/%d', v, N));
    B_all(v,:) = Bv(1:ds:end);
end

%% 8. Create input matrices
Xs  = B_all(:,1:end-1);      % 40 x 303
Xsp = B_all(:,2:end);        % 40 x 303

%% 9. aDMDc
[A, B, F] = aDMDc(Xs,Xsp,stimulus,p);

%% 10. Compare

% ---- (a) Dynamics: estimated A (40x40) vs. ground-truth network ----
C_full = C;
C_full(1:N+1:end) = w_EE;   % put self-excitation on the diagonal, so
                            % C_full plays the same role A's diagonal does

blocks = {'V1->V1 (lateral)',      1:n1,     1:n1;
          'V2->V2 (lateral)',      n1+1:N,   n1+1:N;
          'V1->V2 (feedforward)',  n1+1:N,   1:n1;
          'V2->V1 (feedback)',     1:n1,     n1+1:N};

fprintf('%-22s %10s %10s\n', 'Block', 'mean(A)', 'mean(C)');
for k = 1:size(blocks,1)
    rows = blocks{k,2}; cols = blocks{k,3};
    Ablk = A(rows,cols);
    Cblk = C_full(rows,cols);
    if isequal(rows,cols)
        % same-area block: exclude the diagonal (pure self-terms),
        % keep only lateral coupling between *different* nodes
        Ablk = Ablk(~eye(numel(rows)));
        Cblk = Cblk(~eye(numel(rows)));
    end
    fprintf('%-22s %10.4f %10.4f\n', blocks{k,1}, mean(Ablk(:)), mean(Cblk(:)));
end

figure;
subplot(1,2,1); imagesc(A);      axis square; colorbar; title('Estimated A');
subplot(1,2,2); imagesc(C_full); axis square; colorbar; title('Ground truth (C + w_{EE})');

% ---- (b) Space: estimated pRF (F) vs. true V1 receptive fields ----
n_px = p.n_px;
corr_v1 = nan(n1,1);
for v = 1:n1
    true_img = exp(-((X - p.prf(v).x0).^2 + (Y - p.prf(v).y0).^2) / (2*p.prf(v).sigma^2));
    true_img = true_img / sum(true_img(:));
    est_img  = reshape(F(v,:), n_px, n_px);
    corr_v1(v) = corr(true_img(:), est_img(:));
end

figure;
histogram(corr_v1, 10);
xlabel('correlation (true vs. estimated V1 pRF)'); ylabel('# nodes');
title(sprintf('V1 pRF recovery, median r = %.2f', median(corr_v1)));

fprintf('V1 pRF recovery: median r = %.3f (range %.3f to %.3f)\n', ...
    median(corr_v1), min(corr_v1), max(corr_v1));