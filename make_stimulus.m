function I_ext = make_stimulus(C, t_vec, p)
% MAKE_STIMULUS  One-time external input into V1 only (nA)
%
%   C     : (N x N) connectivity, used only to get the network size
%   t_vec : (1 x n_steps) time vector (s)
%   p     : needs J_ext, mu0, stim_onset, stim_dur
%
%   Returns I_ext : (N x n_steps), nonzero only in row 1 (V1)


N   = size(C, 1);

%index for V1
iV1 = 1;

t_on = t_vec >= p.stim_onset & t_vec < p.stim_onset + p.stim_dur;

%Making the matrix for I_ext (nodes x time steps)
I_ext = zeros(N, numel(t_vec));

% Calculating only in V1 where its applied
I_ext(iV1, t_on) = p.J_ext * p.mu0;
end