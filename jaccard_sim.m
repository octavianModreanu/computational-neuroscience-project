function js = jaccard_sim(est, tru)
% JACCARD_SIM  Jaccard similarity of two weight vectors, as in the
% Connective Field Modelling chapter: sum(min)/sum(max) after clipping
% negative weights to 0 and scaling both vectors to a peak of 1
est = max(est, 0);
tru = max(tru, 0);
if max(est) == 0 || max(tru) == 0
    js = 0;
    return
end
est = est / max(est);
tru = tru / max(tru);
js = sum(min(est, tru)) / sum(max(est, tru));
end
