function u = dmf_input(S, C, w_EE, p, I_ext)
% DMF_INPUT  Computes total input current to each region (Eq. 22b)

%   S : (N x 1) vector of synaptic gating variables
%   C : (N x N) structural connectivity matrix
%   p : parameter struct
%   I_ext : (N x 1) external current in nA

%I_ext given default value of 0 if nothing else passed in
if nargin < 5 || isempty(I_ext)
    I_ext = zeros(size(S));
end

local_term   = w_EE .* S;             % self-excitation
network_term = p.G .* (C * S);          % input from other regions


u = p.J .* (local_term + network_term) + p.I0 + I_ext;

end