function R = run_range_robustness(T, p0, forcing, obs2002, horizons, windowLength, indicatorName)

%
% T             : interval table (from load_parameter_intervals)
% p0            : nominal parameter struct
% forcing       : full multi-year struct (from load_forcing)
% obs2002       : struct from observed_2002()
% horizons      : vector of L (days); use 21 for the main illustration
% windowLength  : W (e.g. 2)
% indicatorName : 'ENEI' (primary); can also run 'ENEI_Tonly','DD'
%
% OUTPUT R: table with the contract columns (continuous + timing endpoints).

vars = {'parameter_id','scenario','parameter_value', ...
        'endpoint_enei_aug1','endpoint_n_peak','nymph_peak_date','biofix_date', ...
        'warning_date_low','warning_date_med','warning_date_high', ...
        'lead_low','lead_med','lead_high', ...
        'lead_error_low','lead_error_med','lead_error_high','status'};
R = cell2table(cell(0,numel(vars)),'VariableNames',vars);

riskLevels = {'low','med','high'};
scnNames = {'low','central','high'};

for i = 1:height(T)
    pid = char(T.parameter_id(i));
    if strcmp(pid,'T_opt'), continue; end
    if lower(string(T.independently_varied(i))) ~= "yes", continue; end
    lo = str2double(string(T.lower(i)));
    hi = str2double(string(T.upper(i)));
    if isnan(lo) || isnan(hi), continue; end
    nomVal = p0.(pid);
    scnVals = [lo, nomVal, hi];

   

    for s = 1:3
        p = apply_parameter_scenario(p0, pid, scnVals(s), forcing.T);   % coupling rules + fT_norm recompute

        % --- 2002 continuous endpoints 
        F2002 = year_forcing(forcing, 2002);
        sim = simulate_enei_model(p, F2002);
        Q = compute_continuous_endpoints(sim);

        % --- indicator trajectory used for WARNING detection
        indTraj = compute_indicator(indicatorName, p, F2002);

        % recalibrate per risk on this scenario (Protocol A) and score vs observed
        wdate = struct('low',NaT,'med',NaT,'high',NaT);
        lead  = struct('low',NaN,'med',NaN,'high',NaN);
        lerr  = struct('low',NaN,'med',NaN,'high',NaN);
        status = "ok";
        for L = horizons  % typically a single L for the main table
            for ir = 1:numel(riskLevels)
                rl = riskLevels{ir};
                [thr, st] = recalibrate_threshold(p, forcing, 2002, rl, L, windowLength, indicatorName);
                if isnat(obs2002.crossings.(rl))
                    status = "no_target_crossing"; continue
                end
                if isnan(thr)
                    status = st; continue
                end
                idx = find(indTraj.value >= thr, 1, 'first');   % score on this indicator
                if isempty(idx)
                    status = "warning_not_reached"; continue
                end
                wdate.(rl) = indTraj.dates(idx);
                lead.(rl)  = days(obs2002.crossings.(rl) - wdate.(rl));
                lerr.(rl)  = lead.(rl) - L;
            end
        end

        R = [R; {pid, scnNames{s}, scnVals(s), Q.ENEI_1Aug, Q.Nmax, ...
                 sim.tNpeak, sim.biofix, wdate.low, wdate.med, wdate.high, ...
                 lead.low, lead.med, lead.high, lerr.low, lerr.med, lerr.high, ...
                 status}]; %#ok<AGROW>
    end
end
end

function [thr, status] = recalibrate_threshold(p, forcing, yy, rl, L, W, indicatorName)
% Protocol A: perturbed parameter used in calibration years too; target
% crossings for calibration years are MODEL-derived (internal), 2002 target
% scored separately by the caller against observed crossings.
status = "ok";
thr = NaN;
calYears = (yy-W):(yy-1);
thetas = nan(1,numel(calYears));
riskVal = risk_value(rl);
for c = 1:numel(calYears)
    Fy = year_forcing(forcing, calYears(c));
    if isempty(Fy.T), thetas(c) = NaN; continue; end
    ind = compute_indicator(indicatorName, p, Fy);
    % model-derived crossing for calibration year
    sim = simulate_enei_model(p, Fy);
    kx = find(sim.N / p.kN_scale >= riskVal, 1, 'first');  % scaled nymph crossing
    if isempty(kx), thetas(c) = NaN; continue; end
    tcross = sim.dates(kx);
    tEval = tcross - days(L);
    if tEval < ind.biofix, thetas(c) = NaN; continue; end
    thetas(c) = interp_value(ind, tEval);
end
valid = thetas(~isnan(thetas));
if isempty(valid)
    status = "insufficient_calibration_years"; return
end
thr = median(valid);
end

function v = risk_value(rl)
switch rl, case 'low', v=0.5; case 'med', v=1.0; case 'high', v=2.0; end
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
