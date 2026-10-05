function [prf, n_voxels] = define_v1_retinotopy(n_x, n_y, extent_x, extent_y, sigma, n_exp)
% DEFINE_V1_RETINOTOPY  Assigns each V1 voxel a receptive field position
%                        on a plain Cartesian grid over the visual field.
%
%   n_x, n_y         : number of grid points along x and y
%                      (total voxels = n_x * n_y)
%   extent_x, extent_y : visual field half-width/half-height (deg) --
%                         grid spans [-extent_x, extent_x] x [-extent_y, extent_y]
%   sigma            : RF width (deg), same for every voxel
%   n_exp            : compressive spatial summation exponent (1 = linear
%                      Gaussian pRF; <1 = compressive, per Kay et al. 2013)
%
%   Returns prf : 1 x n_voxels struct array with fields x0, y0, sigma, n
%                 -- directly usable as p.prf in make_stimulus.m / pRF_response.
%           n_voxels : total voxel count (n_x * n_y)

    x_vals = linspace(-extent_x, extent_x, n_x);
    y_vals = linspace(-extent_y, extent_y, n_y);

    [X, Y] = meshgrid(x_vals, y_vals);
    x0 = X(:);
    y0 = Y(:);

    n_voxels = numel(x0);
    prf(n_voxels) = struct('x0', [], 'y0', [], 'sigma', [], 'n', []);

    for v = 1:n_voxels
        prf(v).x0    = x0(v);
        prf(v).y0    = y0(v);
        prf(v).sigma = sigma;
        prf(v).n     = n_exp;
    end

end