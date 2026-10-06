function [S_t, t, u_t, r_t] = simulate_dmf(C, w_EE, p, T, dt, S0, I_ext, dt_store)
% SIMULATE_DMF  Integrates the DMF model (Gilson et al 2016) with additive 
% white noise directly on dS/dt (Deco et al. 2013 formulation)
%
%   C        : (N x N) structural connectivity matrix
%   p        : parameter struct 
%   T        : total simulation time (s)
%   dt       : integration time step (s) -- e.g., 1e-3
%   S0       : (N x 1) initial condition (optional, default 0.1)
%   I_ext    : external input (nA), optional. Either
%                (N x n_steps)  one column per integration step, or
%                (N x n_frames) one column per stimulus frame (p.frame_dt)
%   dt_store : (optional) outputs are averaged over windows of this length
%              (s); default dt = store every step (original behaviour)
%
%   Returns:
%   S_t : (N x n_store) synaptic gating variable trajectories
%   t   : (1 x n_store) time vector (start of each storage window)
%   u_t : (N x n_store) input current (nA)
%   r_t : (N x n_store) firing rate (Hz)

N = size(C, 1);
n_steps = round(T / dt);

if nargin < 6 || isempty(S0)
    S0 = 0.1 * ones(N, 1);
end
if nargin < 7 || isempty(I_ext)
    I_ext = zeros(N, n_steps);
end
if nargin < 8 || isempty(dt_store)
    dt_store = dt;
end

% storage: one output sample every 'every' integration steps
every   = round(dt_store / dt);
n_store = floor(n_steps / every);
t = (0:n_store-1) * every * dt;

% I_ext given per step, or per frame (look up the frame on screen)
per_step        = size(I_ext, 2) == n_steps;
n_cols          = size(I_ext, 2);
steps_per_frame = round(p.frame_dt / dt);

S_t = zeros(N, n_store);
u_t = zeros(N, n_store);   % input current (nA)
r_t = zeros(N, n_store);   % Firing rate  (hz)
S = S0;
w_EE = w_EE(:);

% running sums over the current storage window
accS = zeros(N, 1);
accU = zeros(N, 1);
accR = zeros(N, 1);

for k = 1:n_steps
    % external input at this step
    if per_step
        Ik = I_ext(:, k);
    else
        Ik = I_ext(:, min(floor((k-1) / steps_per_frame) + 1, n_cols));
    end

    % calculate deterministic without noise
    [dSdt, u, h] = dmf_activity_change(S, C, w_EE, p, Ik);

    % calculate stochastic with noise
    eta = randn(N, 1);
    noise_term = p.sigma .* sqrt(dt) .* eta;   % Euler-Maruyama
    
    % update state 
    S = S + dSdt * dt + noise_term;
    S = min(max(S, 0), 1);   % keep gating variable in [0,1]

    % accumulate, store the window average once the window is full
    accS = accS + S;
    accU = accU + u;
    accR = accR + h;
    if mod(k, every) == 0
        j = k / every;
        S_t(:, j) = accS / every;
        u_t(:, j) = accU / every;
        r_t(:, j) = accR / every;
        accS(:) = 0;
        accU(:) = 0;
        accR(:) = 0;
    end
end

end
