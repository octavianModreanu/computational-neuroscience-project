function [T_scores, T_summary, T_prf] = summarize_scenarios(results, scenarios, conn_names)
% SUMMARIZE_SCENARIOS  Tables of the CF and pRF scores across scenarios.
%
%   T_scores  : one row per scenario x connection (on, CF r, null r, JS)
%   T_summary : per connection, mean scores when ON vs OFF across scenarios
%   T_prf     : per scenario, V1 pRF recovery (median / min / max r)

nS = numel(scenarios);
nC = numel(conn_names);

%% One row per scenario x connection
Scenario = {}; Connection = {}; On = []; CF_r = []; Null_r = []; CF_JS = [];
for s = 1:nS
    for k = 1:nC
        Scenario{end+1,1}   = scenarios{s};          %#ok<AGROW>
        Connection{end+1,1} = conn_names{k};         %#ok<AGROW>
        On(end+1,1)     = results(s).cf(k,1);        %#ok<AGROW>
        CF_r(end+1,1)   = results(s).cf(k,2);        %#ok<AGROW>
        Null_r(end+1,1) = results(s).cf(k,3);        %#ok<AGROW>
        CF_JS(end+1,1)  = results(s).cf(k,4);        %#ok<AGROW>
    end
end
T_scores = table(Scenario, Connection, On, CF_r, Null_r, CF_JS);

%% ON vs OFF per connection
Mean_r_ON   = zeros(nC,1);
Mean_r_OFF  = zeros(nC,1);
Diff_r      = zeros(nC,1);
Mean_JS_ON  = zeros(nC,1);
Mean_JS_OFF = zeros(nC,1);
Scen_ON  = cell(nC,1);
Scen_OFF = cell(nC,1);
for k = 1:nC
    rows = strcmp(T_scores.Connection, conn_names{k});
    on  = rows & T_scores.On == 1;
    off = rows & T_scores.On == 0;
    Mean_r_ON(k)   = mean(T_scores.CF_r(on));
    Mean_r_OFF(k)  = mean(T_scores.CF_r(off));
    Mean_JS_ON(k)  = mean(T_scores.CF_JS(on));
    Mean_JS_OFF(k) = mean(T_scores.CF_JS(off));
    Diff_r(k) = Mean_r_ON(k) - Mean_r_OFF(k);
    Scen_ON{k}  = strjoin(T_scores.Scenario(on)', ',');
    Scen_OFF{k} = strjoin(T_scores.Scenario(off)', ',');
end
Connection = conn_names(:);
T_summary = table(Connection, Scen_ON, Scen_OFF, Mean_r_ON, Mean_r_OFF, Diff_r, ...
                  Mean_JS_ON, Mean_JS_OFF);

%% V1 pRF recovery per scenario
Scenario     = scenarios(:);
PRF_median_r = arrayfun(@(r) median(r.prf_r), results(:));
PRF_min_r    = arrayfun(@(r) min(r.prf_r),    results(:));
PRF_max_r    = arrayfun(@(r) max(r.prf_r),    results(:));
T_prf = table(Scenario, PRF_median_r, PRF_min_r, PRF_max_r);
end
