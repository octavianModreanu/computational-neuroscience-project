function p = dmf_get_params()
% DMF_GET_PARAMS  Returns DMF model parameters (Gilson 2016, Table 3)

p.a        = 270;     % (n/C)
p.b        = 108;     % (Hz)
p.c        = 0.154;   % (s)   (called d in Deco et al. 2013)
p.beta     = 0.641;   % kinetic parameter
p.J        = 0.261;   % scaling constant
p.tau_NMDA = 100e-3;  % NMDA decay time constant (s)
p.I0       = 0.3;     % background input (nA)
p.G        = 1.0;     % global coupling (was 2.5)

% --- connectivity profile (used by make_scenario_connectivity.m) ---
p.g_ff     = 1.0;     % feedforward strength  V1 -> V2
p.g_fb     = 0.5;     % feedback strength     V2 -> V1
p.g_lat    = 0.5;     % lateral (within-area) strength
p.w_EE     = 0.5;     % self-excitation when "lateral" is ON (0 when OFF)

% --- noise ---
p.sigma    = 0.01;    % std of noise on the input current (nA)

% --- external stimulus (Wong & Wang 2006; into V1 only) ---
p.J_ext      = 5.2e-4;   % AMPA coupling of external input (nA/Hz)
p.mu0        = 30;       % stimulus strength (Hz)
end
