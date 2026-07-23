function R = run_range_robustness_B(T, p0, forcing, obs2002, L, W, indicatorName)
% RUN_RANGE_ROBUSTNESS_B  PRIMARY robustness analysis (Protocol B: fixed
% nominal thresholds). 
%
% INPUTS as in run_range_robustness, but L and W are scalars (main setting
% L=21, W=2). indicatorName: 'ENEI' or 'DD'.
%
% OUTPUT R: one row per (parameter, scenario) with fixed-threshold lead times
% and their shift from nominal.

riskLevels = {'low','med','high'};
scnNames   = {'low','central','high'};

% ---- 1. nominal fixed thresholds, calibrated once  ----
yrNom = build_year_results(p0, forcing, obs2002, indicatorName);
thetaFix = struct();
for ir = 1:numel(riskLevels)
    rl = riskLevels{ir};
    [thetaFix.(rl), ~] = calibrate_fixed(yrNom, 2002, rl, L, W);
end

% ---- 2. nominal lead times (reference for the shift) ----
F2002 = year_forcing(forcing, 2002);
indNom = compute_indicator(indicatorName, p0, F2002);
leadNom = struct();
for ir = 1:numel(riskLevels)
    rl = riskLevels{ir};
    leadNom.(rl) = score_lead(indNom, thetaFix.(rl), obs2002.crossings.(rl));
end

% ---- 3. perturb each parameter; thresholds stay fixed ----
vars = {'parameter_id','scenario','parameter_value', ...
        'endpoint_enei_aug1','endpoint_n_peak','nymph_peak_date','biofix_date', ...
        'lead_low','lead_med','lead_high', ...
        'dlead_low','dlead_med','dlead_high','max_abs_dlead','status'};
R = cell2table(cell(0,numel(vars)),'VariableNames',vars);

for i = 1:height(T)
    pid = char(T.parameter_id(i));
    if strcmp(pid,'T_opt'), continue; end
    if lower(string(T.independently_varied(i))) ~= "yes", continue; end
    lo = str2double(string(T.lower(i)));
    hi = str2double(string(T.upper(i)));
    if isnan(lo) || isnan(hi), continue; end
    scnVals = [lo, p0.(pid), hi];

    for s = 1:3
        p = apply_parameter_scenario(p0, pid, scnVals(s), forcing.T);
        sim   = simulate_enei_model(p, F2002);
        Q     = compute_continuous_endpoints(sim);
        indP  = compute_indicator(indicatorName, p, F2002);

        lead = struct(); dlead = struct(); status = "ok";
        for ir = 1:numel(riskLevels)
            rl = riskLevels{ir};
            if isnat(obs2002.crossings.(rl))
                status = "no_target_crossing"; lead.(rl)=NaN; dlead.(rl)=NaN; continue
            end
            lv = score_lead(indP, thetaFix.(rl), obs2002.crossings.(rl));
            lead.(rl) = lv;
            if isnan(lv)
                status = "warning_not_reached"; dlead.(rl)=NaN;
            else
                dlead.(rl) = lv - leadNom.(rl);
            end
        end
        dvals = [dlead.low dlead.med dlead.high];
        maxabs = max(abs(dvals(~isnan(dvals))));
        if isempty(maxabs), maxabs = NaN; end

        R = [R; {pid, scnNames{s}, scnVals(s), Q.ENEI_1Aug, Q.Nmax, ...
                 sim.tNpeak, sim.biofix, lead.low, lead.med, lead.high, ...
                 dlead.low, dlead.med, dlead.high, maxabs, status}]; %#ok<AGROW>
    end
end
end

% ---- helpers ----
function yr = build_year_results(p, forcing, obs2002, indicatorName)
years = 2000:2002;
yr(numel(years)) = struct('year',[],'indicator',[],'crossing',[]);
for k = 1:numel(years)
    y = years(k);
    Fy = year_forcing(forcing, y);
    s  = simulate_enei_model(p, Fy);
    ind = compute_indicator(indicatorName, p, Fy);
    cr = struct();
    for r = ["low","med","high"]
        if y == 2002
            cr.(r) = obs2002.crossings.(r);
        else
            kx = find(s.N / p.kN_scale >= risk_value(r), 1, 'first');
            if isempty(kx), cr.(r) = NaT; else, cr.(r) = s.dates(kx); end
        end
    end
    yr(k) = struct('year',y,'indicator',ind,'crossing',cr);
end
end

function [thr, status] = calibrate_fixed(yr, yy, rl, L, W)
status = "ok"; thr = NaN;
calYears = (yy-W):(yy-1);
thetas = nan(1,numel(calYears));
byYear = containers.Map([yr.year], num2cell(1:numel(yr)));
for c = 1:numel(calYears)
    if ~isKey(byYear, calYears(c)), continue; end
    Y = yr(byYear(calYears(c)));
    tcross = Y.crossing.(rl);
    if isnat(tcross), continue; end
    tEval = tcross - days(L);
    if tEval < Y.indicator.biofix, continue; end
    thetas(c) = interp_value(Y.indicator, tEval);
end
valid = thetas(~isnan(thetas));
if isempty(valid), status = "insufficient_calibration_years"; return; end
thr = median(valid);
end

function lead = score_lead(ind, thr, tcross)
lead = NaN;
if isnan(thr) || isnat(tcross), return; end
idx = find(ind.value >= thr, 1, 'first');
if isempty(idx), return; end
lead = days(tcross - ind.dates(idx));
end

function v = risk_value(rl)
switch char(rl), case 'low', v=0.5; case 'med', v=1.0; case 'high', v=2.0; end
end

function Fy = year_forcing(F, y)
idx = year(F.dates)==y;
Fy = struct('T',F.T(idx),'H',F.H(idx),'R',F.R(idx),'W',F.W(idx), ...
            'dates',F.dates(idx),'year',y);
end

function v = interp_value(ind, t)
idx = find(ind.dates <= t, 1, 'last');
if isempty(idx), v = NaN; return; end
if idx >= numel(ind.dates), v = ind.value(end); return; end
w = days(t-ind.dates(idx))/days(ind.dates(idx+1)-ind.dates(idx));
v = (1-w)*ind.value(idx) + w*ind.value(idx+1);
end
