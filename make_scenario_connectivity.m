function [C, w_EE] = make_scenario_connectivity(scenario, p)
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
    %            ff  fb  lat1 lat2
    case 'A', f = [1   1   0    0];
    case 'B', f = [0   0   1    1];
    case 'C', f = [1   1   0    1];
    case 'D', f = [1   0   1    1];
    case 'E', f = [0   1   1    1];
    case 'F', f = [1   1   1    0];
    otherwise
        error('Unknown scenario "%s". Use A-F.', scenario);
end

% Connectivity mtx is decided by the respective scenario
% the parameters added are feedforward/back strengths of the regions
ff = f(1); fb = f(2); lat1 = f(3); lat2 = f(4);

C = [ 0              fb * p.g_fb ;    % into V1: from V1, from V2
      ff * p.g_ff    0           ];   % into V2: from V1, from V2

%w_EE is the self-excitation weight
w_EE = p.w_EE * [lat1; lat2];        % column vector, one value per nod e

end
