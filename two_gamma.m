function h = two_gamma(t)
% TWO_GAMMA  Canonical two-gamma HRF (same as HGR.two_gamma), t in s
h = (6*t.^5.*exp(-t))./gamma(6) - 1/6*(16*t.^15.*exp(-t))/gamma(16);
end
