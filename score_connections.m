function tab = score_connections(A, W, n1, N)
% SCORE_CONNECTIONS  Compares the rows of the estimated A with the true
% (unscaled) connection weights W, per connection type.
%
%   A  : (N x N) estimated connectivity (aDMDc)
%   W  : weight struct from make_scenario_connectivity
%
%   Returns tab : one row per connection [ff; fb; latV1; latV2] with
%                 [on, median CF r, median null r, median Jaccard].
%   Null = same rows matched to the wrong receiving nodes. Compares SHAPES
%   (A and C are in different units).

blk = {n1+1:N, 1:n1,   W.ff,     W.flags(2,1);
       1:n1,   n1+1:N, W.fb,     W.flags(1,2);
       1:n1,   1:n1,   W.lat_v1, W.flags(1,1);
       n1+1:N, n1+1:N, W.lat_v2, W.flags(2,2)};

tab = zeros(size(blk,1), 4);
for k = 1:size(blk,1)
    rows = blk{k,1};
    cols = blk{k,2};
    Wt   = blk{k,3};

    Ablk = A(rows, cols);
    if isequal(rows, cols)
        Ablk(1:numel(rows)+1:end) = 0;          % ignore self-connections
    end
    Wnull = circshift(Wt, 7, 1);

    r  = zeros(numel(rows), 1);
    nl = r;
    js = r;
    for i = 1:numel(rows)
        r(i)  = corr(Ablk(i,:)', Wt(i,:)');
        nl(i) = corr(Ablk(i,:)', Wnull(i,:)');
        js(i) = jaccard_sim(Ablk(i,:), Wt(i,:));
    end
    tab(k,:) = [blk{k,4}, median(r), median(nl), median(js)];
end
end
