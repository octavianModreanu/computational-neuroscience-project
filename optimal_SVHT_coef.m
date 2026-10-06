function omega = optimal_SVHT_coef(beta)
% OPTIMAL_SVHT_COEF  Gavish & Donoho (2014) optimal hard threshold
% coefficient, noise level unknown:
%   threshold = omega(beta) * median(singular values)
w = (8 * beta) ./ (beta + 1 + sqrt(beta.^2 + 14 * beta + 1));
lambda_star = sqrt(2 * (beta + 1) + w);
omega = lambda_star ./ sqrt(median_marcenko_pastur(beta));
end

function med = median_marcenko_pastur(beta)
% Median of the Marcenko-Pastur distribution (bisection)
lobnd = (1 - sqrt(beta))^2;
hibnd = (1 + sqrt(beta))^2;
topSpec = hibnd;
botSpec = lobnd;
dens   = @(x) sqrt(max((topSpec - x).*(x - botSpec), 0)) ./ (beta .* x) ./ (2*pi);
MarPas = @(x0) 1 - integral(dens, x0, topSpec);
change = 1;
while change && (hibnd - lobnd > .001)
    change = 0;
    xs = linspace(lobnd, hibnd, 5);
    ys = arrayfun(MarPas, xs);
    if any(ys < 0.5), lobnd = max(xs(ys < 0.5)); change = 1; end
    if any(ys > 0.5), hibnd = min(xs(ys > 0.5)); change = 1; end
end
med = (hibnd + lobnd) / 2;
end
