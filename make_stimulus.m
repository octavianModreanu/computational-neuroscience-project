function I_ext = make_stimulus(N, t_vec, p, A, X, Y)
% MAKE_STIMULUS  External input into V1, driven by a bar stimulus seen
% through V1's population receptive field (pRF).
%
%   N     : number of nodes
%   t_vec : (1 x n_steps) simulation time vector (s)
%   p     : parameter struct; needs J_ext, mu0, frame_dt, and prf(1)
%           with fields x0, y0, sigma (deg) and n (exponent, 1 = linear)
%   A     : (n_pixels x n_frames) stimulus, one flattened frame per column
%   X, Y  : (150 x 150) visual-field coordinates of each pixel (deg)
%
%   Returns I_ext : (N x n_steps) input current (nA), nonzero only in V1


n_frames = size(A, 2);

%index for V1
iV1 = 1;


% Overlap between the bar and V1's pRF on every frame, in [0, 1]:
% 0 = bar misses the pRF, 1 = bar covers all of it
r = pRF_response(A, X, Y, p.prf(iV1));             % 1 x n_frames


% Which stimulus frame is on screen at each simulation step.
% Each frame lasts frame_dt seconds (one TR), i.e. many integration steps.
idx = min(floor(t_vec / p.frame_dt) + 1, n_frames);
assert(idx(1) == 1 && idx(end) == n_frames, ...
    'Simulation length does not match the stimulus length.');


% Convert overlap into current: full coverage gives J_ext * mu0 nA
I_ext = zeros(N, numel(t_vec));

I_ext(iV1, :) = p.J_ext * p.mu0 * r(idx);
end