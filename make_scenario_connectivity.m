function [C, w_EE, W] = make_scenario_connectivity(scenario, p, n1, n2, N, anchor)
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
%   W        : unscaled weight matrices (ground-truth CFs, used for scoring)
%              W.ff (n2 x n1), W.fb (n1 x n2), W.lat_v1 (n1 x n1),
%              W.lat_v2 (n2 x n2), W.flags (= Mflag), W.prf_pos_v2 (n2 x 2)
%
%   Scenario   feedforward  feedback  lateral V1  lateral V2
%      A           x           x
%      B                                  x           x
%      C           x           x                      x
%      D           x                      x           x
%      E                       x          x           x
%      F           x           x          x

prf_pos_v1 = [[p.prf.x0]',[p.prf.y0]'];


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

%Rows are receivers, cols are senders

%% Feedforward connections
% Each v2 node is made up by a sampled portion of the v1 nodes. Each v2
% samples an anchor at random from v1, then a gaussian around the anchor
% connects to the v2 node, making up its prf.


% Each row = one v2 node;
% cols: x,y coordinates of the v2 node's anchor in v1

v2_anchor_pos = prf_pos_v1(anchor, :);

% distance between each anchor to every v1 node

dx_ff = v2_anchor_pos(:,1) - prf_pos_v1(:,1)';
dy_ff = v2_anchor_pos(:,2) - prf_pos_v1(:,2)';
dist_ff = dx_ff.^2 + dy_ff.^2;

% gaussian falloff;
% Each row = one v2 node
% Each col = one v1 node
% Each entry W(i,j) = how strongly connected each v1 node (j) is connected
% to the (i) v2 node

W_ff = exp(-dist_ff/(2 * p.s_ff^2));
W_ff = W_ff ./ sum(W_ff, 2);

prf_pos_v2 = W_ff * prf_pos_v1; % (much larger) v2 prfs

%% Feedback connections
% The feedback connection is based on how close the position of each v1 node 
% is to v2 positions. If the prf of a node in v2 and v1 are the same, then
% they should be strongly connected. 
% 
% Notes:
% ( idk how biological this is tbh, it'd also make sense for feedback connections
% to be sent back to the feedforward positions but ig this is what's
% already doing :shrug: )
%
% 


dx_fb = prf_pos_v1(:,1) - prf_pos_v2(:,1)';
dy_fb = prf_pos_v1(:,2) - prf_pos_v2(:,2)';
dist_fb = dx_fb.^2 + dy_fb.^2;

% Each v1 node i gets its feedback input from a v2 node whose position on the
% visual field is close to v1 node i's own visual field position
W_fb = exp(-dist_fb/(2*p.s_fb^2));


%% Lateral connections

dx_lat_v1 = prf_pos_v1(:,1) - prf_pos_v1(:,1)';
dy_lat_v1 = prf_pos_v1(:,2) - prf_pos_v1(:,2)';
dist_lat_v1 = dx_lat_v1.^2 + dy_lat_v1.^2;

W_lat_v1 = exp(-dist_lat_v1/(2*p.s_lat^2));

dx_lat_v2 = prf_pos_v2(:,1) - prf_pos_v2(:,1)';
dy_lat_v2 = prf_pos_v2(:,2) - prf_pos_v2(:,2)';
dist_lat_v2 = dx_lat_v2.^2 + dy_lat_v2.^2;

W_lat_v2 = exp(-dist_lat_v2/(2*p.s_lat^2));


%% Normalize and build C

nrm = @(W) W ./ sum(W, 2);

W_fb = nrm(W_fb);

W_lat_v1(1:n1+1:end) = 0;
W_lat_v1 = nrm(W_lat_v1);

W_lat_v2(1:n2+1:end) = 0;
W_lat_v2 = nrm(W_lat_v2);

i1 = 1:n1;  i2 = n1+1:N;
C = zeros(N);
C(i1,i1) = Mflag(1,1) * Gain(1,1) * W_lat_v1;
C(i1,i2) = Mflag(1,2) * Gain(1,2) * W_fb;
C(i2,i1) = Mflag(2,1) * Gain(2,1) * W_ff;
C(i2,i2) = Mflag(2,2) * Gain(2,2) * W_lat_v2;

w_EE = p.w_EE * ones(N,1);

%% Unscaled weights and flags (ground truth for scoring)
W.ff         = W_ff;
W.fb         = W_fb;
W.lat_v1     = W_lat_v1;
W.lat_v2     = W_lat_v2;
W.flags      = Mflag;
W.prf_pos_v2 = prf_pos_v2;
