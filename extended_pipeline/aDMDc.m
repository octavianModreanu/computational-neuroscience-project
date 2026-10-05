function [A, B, F] = aDMDc(X, Xp, stimulus, p)
% aDMDc - algebraic dynamic mode decomposition with control
% this function takes in X and Xp, as well as the stimulus, and hashes the
% stimulus. then it performs aDMDc. 
% this function returns A, B, and F, the matrices defined by equations
% 3.19-3.21 of the chapter from the PhD we were given. 
% 
arguments (Input)
    X
    Xp
    stimulus 
    p
end

arguments (Output)
    A
    B
    F
end

addpath("CNI_toolbox/code/matlab/")    % for HGR, same as retinotopy.m
parameters.fwhm         = 0.15;
parameters.eta          = 0.1;
parameters.f_sampling   = 1/p.frame_dt;   % 1/TR
parameters.r_stimulus   = p.n_px;         % 150
parameters.n_features   = 250;
parameters.n_gaussians  = 5;
parameters.n_voxels     = 2;              % V1 and V2
prf_mapper = HGR(parameters);
gam = prf_mapper.get_features;            % 22500 x 250

shash = (stimulus' * gam)';               % 250 x 304 (hashed stimulus, pixel space -> feature space)
Ups  = zscore(shash(:,1:end-1)')';
dUps = zscore(shash(:,2:end)')';

Xz  = zscore(X')';
Xpz = zscore(Xp')';
Omega = [Xz; Ups; dUps];

[U,Sig,V] = svd(Omega, 'econ');
beta = size(Sig,2)/size(Sig,1);
threshold = optimal_SVHT_coef(beta,Sig);
rtil = sum(diag(Sig) >= threshold);
Util = U(:,1:rtil); Sigtil = Sig(1:rtil,1:rtil); Vtil = V(:,1:rtil);
disp(rtil)

n = size(X,1); q = size(Ups,1);
U_1 = Util(1:n,:);
U_2 = Util(n+1:n+q,:);
U_3 = Util(n+q+1:n+q+q,:);

A = Xpz * Vtil * inv(Sigtil) * U_1';   % 
B = Xpz * Vtil * inv(Sigtil) * U_2';   %
F = Xpz * Vtil * inv(Sigtil) * U_3';   %

F = (gam * F')';    % 2 x 22500 — NOW back in pixel space, one row per region
end