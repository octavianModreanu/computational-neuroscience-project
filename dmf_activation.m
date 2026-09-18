function h = dmf_activation(x, p)
% DMF_ACTIVATION  Sigmoidal input-output function (Eq. 22, H(x))
%
%   x : input current (vector, one value per region)
%   p : parameter struct (needs a, b, d)

num = p.a .* x - p.b;
den = 1 - exp(-p.c .* (p.a .* x - p.b));

h = num ./ den;

% Handle the 0/0 singularity as x -> b/a (L'Hopital limit = 1/d)
h(den == 0) = 1 / p.c;

end