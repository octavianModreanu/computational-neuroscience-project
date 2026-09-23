function [C, w_EE] = make_scenario_connectivity(scenario,area, p)
% MAKE_SCENARIO_CONNECTIVITY  Ground-truth connectivity for scenarios A-F
%
%   Two nodes: node 1 = V1 (lower area), node 2 = V2 (higher area)
%
%   Convention: C(i, j) = connection FROM node j INTO node i
%       C(2,1) = feedforward (V1 -> V2)
%       C(1,2) = feedback    (V2 -> V1)
%   Lateral connections = self-excitation, stored in w_EE (N x 1)
%
%   scenario : 'A' ... 'F'
%   p        : parameter struct (needs g_ff, g_fb, w_lat)
%
%   Scenario   feedforward  feedback  lateral V1  lateral V2
%      A           x           x
%      B                                  x           x
%      C           x           x                      x
%      D           x                      x           x
%      E                       x          x           x
%      F           x           x          x


% in this part, each case makes a vector f where ff is V1
% fb is V2, lat1 is v1 self-excitation, lat2 is v2 self-excitation
switch upper(scenario)
    case 'A', Mflag = [0 1; 1 0];
    case 'B', Mflag = [1 0; 0 1];
    case 'C', Mflag = [0 1; 1 1];
    case 'D', Mflag = [1 0; 1 1];
    case 'E', Mflag = [1 1; 0 1];
    case 'F', Mflag = [1 1; 1 0];
    otherwise, error('Unknown scenario "%s". Use A-F.', scenario);
end
% Connectivity mtx is decided by the respective scenario
% the parameters added are feedforward/back strengths of the regions
Gain = [p.g_lat  p.g_fb ;
        p.g_ff   p.g_lat];

N = numel(area);
C = Mflag(area, area) .* Gain(area, area);   % expand 2x2 to N x N
C(1:N+1:end) = 0;

% Normalise each block so total input does not grow with N
% C*S in dmf_input sums over all incoming connections, so without this the
% total input to a node would grow with the number of nodes in the sending
% area: at N=40 a V2 node would receive ~20x the feedforward drive it gets
% at N=2, pushing u well past the sigmoid's linear range and pinning every
% node at S=1.
%
% Dividing by the number of senders makes each block an average over its
% source population rather than a sum
for a = 1:2
    for b = 1:2
        ii = area == a;  jj = area == b;
        n  = nnz(jj) - (a == b);             % senders, minus self
        if n > 0, C(ii, jj) = C(ii, jj) / n; end
    end
end

%w_EE is the self-excitation weight
w_EE = p.w_EE * ones(N, 1);

end
