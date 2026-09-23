function I_ext = make_stimulus(N, t_vec, p, A, X, Y)
% MAKE_STIMULUS  External input into V1, driven by a bar stimulus seen
% through each V1 voxel's population receptive field (pRF).
%
%   N     : number of nodes
%   t_vec : (1 x n_steps) simulation time vector (s)
%   p     : parameter struct; needs J_ext, mu0, frame_dt, and prf, a
%           1 x n_v1 struct array (one entry per V1 voxel) with fields
%           x0, y0, sigma (deg) and n (exponent, 1 = linear) -- as
%           returned by define_v1_retinotopy.m
%   A     : (n_pixels x n_frames) stimulus, one flattened frame per column
%   X, Y  : (150 x 150) visual-field coordinates of each pixel (deg)
%
%   Returns I_ext : (N x n_steps) input current (nA), nonzero only in V1
%                   (the first n_v1 = numel(p.prf) rows -- V1 voxels are
%                   assumed to occupy the first block of node indices,
%                   matching the area convention used elsewhere)


n_frames = size(A, 2);

% V1 voxels occupy the first n_v1 node indices
n_v1 = numel(p.prf);


% Which stimulus frame is on screen at each simulation step.
% Each frame lasts frame_dt seconds (one TR), i.e. many integration steps.
% Same for every voxel, so compute once outside the loop.
idx = min(floor(t_vec / p.frame_dt) + 1, n_frames);
assert(idx(1) == 1 && idx(end) == n_frames, ...
    'Simulation length does not match the stimulus length.');


% Convert overlap into current: full coverage gives J_ext * mu0 nA
I_ext = zeros(N, numel(t_vec));

for v = 1:n_v1
    % Overlap between the bar and this voxel's pRF on every frame, in
    % [0, 1]: 0 = bar misses the pRF, 1 = bar covers all of it
    r = pRF_response(A, X, Y, p.prf(v));        % 1 x n_frames

    I_ext(v, :) = p.J_ext * p.mu0 * r(idx);
end

end