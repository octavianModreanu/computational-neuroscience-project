function gam = hgr_features(r_stimulus, n_features, n_gaussians, fwhm)
% HGR_FEATURES  Hashed Gaussian features, identical to HGR.create_gamma
% (CNI toolbox, Bhat et al. 2021): each feature = sum of n_gaussians
% Gaussians at random pixel positions, normalised to sum 1.
%
%   r_stimulus  : image side length (pixels)
%   n_features  : number of features
%   n_gaussians : Gaussians per feature
%   fwhm        : Gaussian FWHM as a fraction of the image side
%
%   Returns gam : (n_pixels x n_features)

Yg = linspace(0, r_stimulus, r_stimulus)' * ones(1, r_stimulus);
Xg = ones(r_stimulus, 1) * linspace(0, r_stimulus, r_stimulus);
sigma = fwhm * r_stimulus / (2 * sqrt(2 * log(2)));

pix_id = linspace(0, r_stimulus^2 - 1, n_features * n_gaussians);
pix_id = pix_id(randperm(n_features * n_gaussians));
x = floor(pix_id / r_stimulus) + 1;
y = mod(pix_id, r_stimulus) + 1;

gam = zeros(r_stimulus^2, n_features);
for i = 0:n_features-1
    g = zeros(r_stimulus);
    for j = 1:n_gaussians
        k = i*n_gaussians + j;
        g = g + exp(-((Xg - x(k)).^2 + (Yg - y(k)).^2) ./ (2 * sigma^2));
    end
    gam(:, i+1) = g(:) / sum(g(:));
end
end
