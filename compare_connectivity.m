function metrics = compare_connectivity(C, A, scenario_label)
% COMPARE_CONNECTIVITY  Compare aDMDc-estimated connectivity A against
% ground-truth structural connectivity C, for a 2-node (V1/V2) system.
%
%   metrics = compare_connectivity(C, A, 'A')
%
% IMPORTANT: C by construction has zero diagonal -- lateral/self
% excitation lives in w_EE, not in C (see make_scenario_connectivity.m).
% A's diagonal instead reflects whatever local dynamics (self-excitation
% AND the -S/tau_NMDA decay term) look like from the outside. So the two
% diagonals are NOT estimating the same physical quantity and should not
% be compared directly. This function reports them separately: off-diagonal
% entries (the actual feedforward/feedback connectivity you care about)
% get compared to C; diagonal entries are reported for inspection only.

    if nargin < 3, scenario_label = ''; end

    ff_true = C(2,1);  fb_true = C(1,2);   % ground truth V1->V2, V2->V1
    ff_est  = A(2,1);  fb_est  = A(1,2);   % estimated equivalents

    metrics.scenario  = scenario_label;
    metrics.ff_true   = ff_true;   metrics.ff_est = ff_est;
    metrics.fb_true   = fb_true;   metrics.fb_est = fb_est;
    metrics.diag_A    = diag(A)';  % reported, not scored against anything

    % Relative error where a connection is actually present; sign check
    % where it's meant to be absent (should stay near 0 / same sign as noise floor)
    metrics.ff_rel_err = safe_rel_err(ff_true, ff_est);
    metrics.fb_rel_err = safe_rel_err(fb_true, fb_est);
    metrics.sign_match = (sign(ff_true) == sign(ff_est)) + (sign(fb_true) == sign(fb_est));
    metrics.sign_match = metrics.sign_match / 2;   % fraction in [0, 0.5, 1]

    fprintf('Scenario %-2s | ff: true=%.3f est=%.3f (rel.err=%.2f) | fb: true=%.3f est=%.3f (rel.err=%.2f) | sign match=%.0f%%\n', ...
        scenario_label, ff_true, ff_est, metrics.ff_rel_err, fb_true, fb_est, metrics.fb_rel_err, 100*metrics.sign_match);
end

function e = safe_rel_err(truth, est)
    if truth == 0
        e = abs(est);         % "relative" error undefined at 0; report the raw leak instead
    else
        e = abs(est - truth) / abs(truth);
    end
end
