function [A, B, F, info] = aDMDc(X, Xp, stimulus, p)
% ADMDC  Algebraic DMD with control:
%           X(:,k+1) = A X(:,k) + B u(:,k) + F u(:,k+1)
%   u = stimulus hashed into HGR features, delayed to match the BOLD lag.
%
%   X, Xp    : (n x m) BOLD snapshots and the same shifted by one frame
%   stimulus : (n_pixels x n_frames) the frames that were simulated; the
%              last m+1 frames are aligned with the BOLD
%   p        : parameter struct (needs n_px, frame_dt, stim_lag, dmd_rank)
%
%   Returns A (n x n), B (n x q), F (n x n_pixels) in original units, and
%   info (rank, svht_rank, threshold, beta, singular values)

%% Hashed stimulus features
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

shash  = shash';                                  % 250 x n_frames
n_snap = size(X,2) + 1;
shash  = shash(:, end-n_snap+1:end);              % same frames as the BOLD
U  = shash(:, 1:end-1);
dU = shash(:, 2:end);

%% SVD of the stacked data (z-scored rows, for the SVD only)
[Xz,  ~, sX ] = zscore(X,  0, 2);
[Xpz, ~, sXp] = zscore(Xp, 0, 2);
[Uz,  ~, sU ] = zscore(U,  0, 2);
[dUz, ~, sdU] = zscore(dU, 0, 2);
sU(sU == 0)   = 1;
sdU(sdU == 0) = 1;

Omega = [Xz; Uz; dUz];
[Uo, Sig, V] = svd(Omega, 'econ');
sv = diag(Sig);

%% Rank cut-off: SVHT (beta from Omega's shape) or forced rank
beta      = min(size(Omega)) / max(size(Omega));
threshold = optimal_SVHT_coef(beta) * median(sv);
svht_rank = sum(sv >= threshold);
if isfield(p, 'dmd_rank') && ~isempty(p.dmd_rank)
    rtil = p.dmd_rank;
else
    rtil = svht_rank;
end
fprintf('aDMDc: Omega %dx%d, beta = %.3f, SVHT rank = %d, used rank = %d\n', ...
    size(Omega,1), size(Omega,2), beta, svht_rank, rtil);

%% Operators
Util   = Uo(:, 1:rtil);
Sigtil = Sig(1:rtil, 1:rtil);
Vtil   = V(:, 1:rtil);
n = size(X, 1);
q = size(U, 1);

K  = Xpz * Vtil / Sigtil;
Az = K * Util(1:n, :)';
Bz = K * Util(n+1:n+q, :)';
Fz = K * Util(n+q+1:n+2*q, :)';

% back to original units
A = sXp .* Az ./ sX';
B = sXp .* Bz ./ sU';
F = sXp .* Fz ./ sdU';
F = (gam * F')';                                  % n x n_pixels

info.rank      = rtil;
info.svht_rank = svht_rank;
info.threshold = threshold;
info.beta      = beta;
info.sv        = sv;
end
