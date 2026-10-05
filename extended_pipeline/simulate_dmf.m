function [S_t, t, u_t, r_t] = simulate_dmf(C, w_EE, p, T, dt, S0, I_ext)
% SIMULATE_DMF  Integrates the DMF model (Gilson et al 2016) with additive 
% white noise directly on dS/dt (Deco et al. 2013 formulation)
%
%   C  : (N x N) structural connectivity matrix
%   p  : parameter struct 
%   T  : total simulation time (s)
%   dt : integration time step (s) -- e.g., 1e-3
%   S0 : (N x 1) initial condition (optional, default 0.1)
%
%   Returns:
%   S_t : (N x n_steps) synaptic gating variable trajectories
%   t   : (1 x n_steps) time vector

N = size(C, 1);
n_steps = round(T / dt);
t = (0:n_steps-1) * dt;

if nargin < 6 || isempty(S0)
    S0 = 0.1 * ones(N, 1);
end
if nargin < 7 || isempty(I_ext)
    I_ext = zeros(N, n_steps);
end

S_t = zeros(N, n_steps);
u_t = zeros(N, n_steps);   % input current (nA)
r_t = zeros(N, n_steps);   % Firing rate  (hz)
S = S0;



for k = 1:n_steps
    % calculate deterministic without noise
    [dSdt,u,h] = dmf_activity_change(S, C, w_EE(:), p, I_ext(:, k));

    u_t(:, k) = u;
    r_t(:, k) = h;

    
    % calculate stochastic with noise
    eta = randn(N, 1);
    noise_term = p.sigma .* sqrt(dt) .* eta;   % Euler-Maruyama
    
    % update state 
    S = S + dSdt * dt + noise_term;
    S = min(max(S, 0), 1);   % keep gating variable in [0,1]

    S_t(:, k) = S;
end

end