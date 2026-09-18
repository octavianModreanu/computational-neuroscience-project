function [dS_dt, u, h] = dmf_activity_change(S, C, w_EE, p, I_ext)
% DMF_RHS  Computes dS/dt for the full network (Eq. 22a)
%
%   S : (N x 1) current state (synaptic gating variables)
%   C : (N x N) structural connectivity matrix
%   p : parameter struct

%if no argument for I_ext is passed, make it empty matrix
if nargin < 5, I_ext = []; 
end

u = dmf_input(S, C,w_EE, p, I_ext);   % calculate the total input
h = dmf_activation(u, p); % calculate the sigmoid function of u

decay = -S ./ p.tau_NMDA;
drive = p.beta .* (1 - S) .* h;

dS_dt = decay + drive;

end