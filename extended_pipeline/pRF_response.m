function r = pRF_response(A, X, Y, prf)
% PRF_RESPONSE  Overlap between each stimulus frame and a Gaussian pRF.
%   Returns r : (1 x n_frames), in [0, 1]



% x0,y0 and sigma are already pre-assigned in main_exec.m --> maybe not
% mistake?
%--------------------------------------------------------------------------
% prf: struct with x0, y0, sigma (deg), n (compressive exponent, 1 = linear)
%making example pRF struct:
% prf.x0    = 3;      % pRF centre, degrees right of fixation
% prf.y0    = -2;     % pRF centre, degrees below fixation
% prf.sigma = 1.5;    % pRF size, degrees

%when we scale up this will obviously change, x0,y0 and sigma will be assigned to
%each node according to something?
%--------------------------------------------------------------------------


%gaussian (eq.2 dumoulin and wandell)
g = exp(-((X - prf.x0).^2 + (Y - prf.y0).^2) / (2*prf.sigma^2));

%normalization (flatten matrix then divide every entry by sum - useful if 
% we end up using different sigmas for different nodes)
%flattening of matrix is also for computing the weighted sum 
g = g(:) / sum(g(:));

% dot product of g'*A, gives a weighted sum, how much of the pRF the bar
%is covering
r = (g' * double(A));      % 1 x n_frames, in [0, 1]
end
