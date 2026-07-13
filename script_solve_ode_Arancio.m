% ==============================================================
% Jacobiasca lybica Stage-Structured Model (per leaf)
% ODE simulation + ENEI early-warning calibration
%
% Thermal case:
%   Analytis:  f_T = (T-Tmin)*sqrt(Tmax-T)
%
% Main features:
% - Feudo Arancio 2002 comparison using peak-based kA and kN
%   (scales the MODEL to the units of the observations)
% - Three-level nymph-risk classification:
%       Low    = 0.5 nymphs/leaf
%       Medium = 1.0 nymphs/leaf
%       High   = 2.0 nymphs/leaf
% - ENEI = cumulative egg-to-nymph emergence pressure
% - ENEI warning thresholds calibrated for L = 7, 14, 21 days
% - Moving-window operational calibration with W = 2 years
% - 2002 warning lead times evaluated against OBSERVED risk-crossing
%   dates obtained by linear interpolation of the bi-weekly field data
% - Degree-day benchmark included for comparison with ENEI
%
% Required input file:
%   Arancio_dati_2000_2025.mat
%   containing daily Tvec, Hvec, Rvec, Wvec from 2000 to 2025
%
% Main outputs:
%   figures_paper/   figures used in the paper
%   tables_paper/    csv tables used in the paper
% ==============================================================

clc; clear; close all

%% 1. USER SETTINGS
% Define input file, target year, warning horizons, risk thresholds and output folders.
dataFile  = 'Arancio_dati_2000_2025.mat';
if ~isfile(dataFile) && isfile('Arancio_dati_2000_2025(2).mat')
    dataFile = 'Arancio_dati_2000_2025(2).mat';
end
anno_cmp  = 2002;
years_eval = 2000:2025;
lead_days  = [7 14 21];
selectedLeadForFigures = 21;

% Nymph-density thresholds [nymphs / leaf]
% Revised three-level operational classification:
%   low    = 0.5 nymphs/leaf
%   medium = 1.0 nymphs/leaf
%   high   = 2.0 nymphs/leaf
N_thr_low  = 0.5;
N_thr_med  = 1.0;
N_thr_high = 2.0;
riskNames = ["Low","Medium","High"];
riskVals  = [N_thr_low, N_thr_med, N_thr_high];

% Rolling hindcast aggregation
thresholdStat = "median";   % options: "median", "q25", "q10"
minPrevYears  = 1;

% Moving window length used for Figures 10-11 (and recommended operationally)
W_movwin = 2;

% Output
saveFigures = true;
saveTables  = true;
figFolder   = 'figures_paper';
tabFolder   = 'tables_paper';
if saveFigures && ~exist(figFolder,'dir'), mkdir(figFolder); end
if saveTables  && ~exist(tabFolder,'dir'), mkdir(tabFolder); end

%% 2. LOAD DAILY METEOROLOGICAL DATA
% The .mat file must contain Tvec, Hvec, Rvec and Wvec with the same length.
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end
load(dataFile)

requiredVars = {'Tvec','Hvec','Rvec','Wvec'};
for iv = 1:numel(requiredVars)
    if ~exist(requiredVars{iv}, 'var')
        error('The data file must contain variable: %s', requiredVars{iv});
    end
end

% Published archive is cleaned; retain a hard guard against unprocessed NASA
% POWER fill values or non-finite entries in replacement datasets.
if any(Tvec(:) <= -900) || any(Hvec(:) <= -900) || ...
        any(Rvec(:) <= -900) || any(Wvec(:) <= -900)
    error('Forcing contains an unprocessed fill value (<= -900).');
end
if any(~isfinite(Tvec(:))) || any(~isfinite(Hvec(:))) || ...
        any(~isfinite(Rvec(:))) || any(~isfinite(Wvec(:)))
    error('Forcing contains NaN or Inf values.');
end

%% 3. TIME GRID
% Daily grid from 1 January 2000. The ODE is advanced with one-day Euler steps.
nt = length(Tvec);
if any([length(Hvec), length(Rvec), length(Wvec)] ~= nt)
    error('Tvec, Hvec, Rvec and Wvec must have the same length.');
end
startDate = datetime(2000,1,1);
dates = startDate + days(0:nt-1);
ht = 1.0;

%% Model parameters -- SINGLE SOURCE AND SHARED RIGHT-HAND SIDE
% All parameter values are loaded from the authoritative CSV and assembled in
% the native parameter structure used by enei_rhs.m. The main simulation,
% time-resolved sensitivities, endpoint elasticities and range analyses now
% evaluate exactly the same population-dynamics routine.
ParamTable = load_parameter_intervals('ENEI_parameter_intervals_updated.csv');
p_model = build_nominal_parameters(ParamTable);
p_model = set_fT_normalization(p_model, Tvec);

% Aliases retained only for downstream calibration/plotting functions.
% They are all read from p_model; no numerical parameter is duplicated here.
beta_max   = p_model.beta_max;
muE_min    = p_model.mu_E_min;
alphaT_E   = p_model.alpha_T;
alphaH_E   = p_model.alpha_H;
gammaE_max = p_model.gamma_E_max;
muN_min    = p_model.mu_N_min;
kW_mort    = p_model.k_W;
gammaN_max = p_model.gamma_N_max;
alphaT_N   = p_model.delta_T;
alphaH_N   = p_model.delta_H;
muA_min    = p_model.mu_A_min;
betaT_A    = p_model.beta_T;
betaH_A    = p_model.beta_H;

%% Environmental factors (from the same single source)
Tmin       = p_model.T_min;
Tmax       = p_model.T_max;
Tmin_spost = p_model.T_move;
Topt       = p_model.T_opt;      % DERIVED: (2*Tmax+Tmin)/3
cH   = p_model.c_H;
H0   = p_model.H_m;
Hopt = p_model.H_m;
cR   = p_model.c_R;
R0   = p_model.R_m;
KW_half    = p_model.K_W;
alphaW_dev = p_model.alpha_W;

fH = @(H) 1.0 ./ (1.0 + exp(-cH*(H - H0)));
fR = @(R) 1.0 ./ (1.0 + exp(-cR*(R - R0)));
fW = @(W) W ./ (W + KW_half);
gW = @(W) 1 ./ (1 + alphaW_dev*W);

%% 6. OBSERVED FIELD DATA FOR 2002
% Adults are trap counts; nymphs are community-level nymphs per leaf.
obsDatesA = datetime(2002, [5 6 6 7 7 7 8 8 9 9 10 10], ...
                           [22 5 19 3 17 31 14 28 11 25 9 23])';
obsAdultsTrap = [0 0 0 10 35 80 150 140 200 265 230 190]';

obsDatesN = datetime(2002, [5 5 6 6 6 7 7 8 8 9 9 10 10], ...
                           [5 19 2 16 30 14 28 11 25 8 22 6 20])';
obsNymphsAllLeaf = [0.00 0.00 0.08 0.12 0.20 0.25 0.27 0.28 0.61 2.00 4.10 2.00 0.43]';

% Interpolated observed risk-crossing dates for 2002.
% These dates are used to evaluate warning lead times against field data
% rather than against the modelled Nscaled trajectory.
obs_tr_2002 = NaT(1, numel(riskVals));
for ir = 1:numel(riskVals)
    thr = riskVals(ir);
    for k = 1:numel(obsNymphsAllLeaf)-1
        if obsNymphsAllLeaf(k) < thr && obsNymphsAllLeaf(k+1) >= thr
            frac = (thr - obsNymphsAllLeaf(k)) / ...
                   (obsNymphsAllLeaf(k+1) - obsNymphsAllLeaf(k));
            delta_days = round(frac * days(obsDatesN(k+1) - obsDatesN(k)));
            obs_tr_2002(ir) = obsDatesN(k) + days(delta_days);
            break
        end
    end
end

fprintf('\nObserved nymph risk crossings (2002, interpolated):\n');
for ir = 1:numel(riskNames)
    if isnat(obs_tr_2002(ir))
        fprintf('  %-7s (%.1f nymphs/leaf): not reached\n', ...
            char(riskNames(ir)), riskVals(ir));
    else
        fprintf('  %-7s (%.1f nymphs/leaf): %s\n', ...
            char(riskNames(ir)), riskVals(ir), datestr(obs_tr_2002(ir),'dd-mmm-yyyy'));
    end
end

%% 7. PLOT SETTINGS
fsAxes = 12;
fsTitle = 14;
fsLabel = 13;
fsLegend = 11;

%% 8. THERMAL RESPONSE FORMULATION
% The raw thermal function is normalized by its maximum value before use.
cases(1).name  = "Analytis";
cases(1).fTraw = @(T) fT_Analytis_raw(T, Tmin, Tmax);
cases(1).Topt  = (2*Tmax + Tmin)/3;

%% 9. TABLE STORAGE
% Tables are progressively filled inside the main loop.
AllValidation = table();
AllThetaYear  = table();
AllRolling    = table();
AllCurves     = table();

%% 10. METEOROLOGICAL FORCING FIGURE FOR 2002
maskFig3 = year(dates) == anno_cmp;
datesF3  = dates(maskFig3);

fig = figure('Color','w','Position',[100 100 1100 640]);
tl  = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile
plot(datesF3, Tvec(maskFig3), 'Color',[0.00 0.45 0.74], 'LineWidth',1.2)
grid on; title('Temperature','FontSize',fsTitle)
xlabel('Date','FontSize',fsLabel); ylabel('T (\circC)','FontSize',fsLabel)
set(gca,'FontSize',fsAxes)

nexttile
plot(datesF3, Hvec(maskFig3), 'Color',[0.85 0.16 0.16], 'LineWidth',1.2)
grid on; title('Relative Humidity','FontSize',fsTitle)
xlabel('Date','FontSize',fsLabel); ylabel('H (%)','FontSize',fsLabel)
set(gca,'FontSize',fsAxes)

nexttile
plot(datesF3, Wvec(maskFig3), 'Color',[0.20 0.70 0.20], 'LineWidth',1.2)
grid on; title('Wind Speed','FontSize',fsTitle)
xlabel('Date','FontSize',fsLabel); ylabel('W (m s^{-1})','FontSize',fsLabel)
set(gca,'FontSize',fsAxes)

nexttile
plot(datesF3, Rvec(maskFig3), 'Color',[0.85 0.20 0.55], 'LineWidth',1.2)
grid on; title('Solar Radiation','FontSize',fsTitle)
xlabel('Date','FontSize',fsLabel); ylabel('R (MJ m^{-2} d^{-1})','FontSize',fsLabel)
set(gca,'FontSize',fsAxes)

title(tl, sprintf('Daily meteorological forcing, Feudo Arancio %d', anno_cmp), ...
    'FontSize',fsTitle)
if saveFigures
    disable_axes_toolbar(fig);
    exportgraphics(fig, fullfile(figFolder,'fig03_meteo_forcing_2002.png'),'Resolution',300);
    exportgraphics(fig, fullfile(figFolder,'fig03_meteo_forcing_2002.pdf'));
end

%% 11. MAIN LOOP OVER THERMAL FORMULATIONS
% For each thermal case: simulate ODE, validate 2002 outputs, calibrate ENEI.
for c = 1:numel(cases)

    thermoName = cases(c).name;
    Topt_case  = cases(c).Topt;

    fprintf('\n====================================================\n');
    fprintf('Thermal formulation: %s\n', char(thermoName));
    fprintf('====================================================\n');

    ftval = zeros(nt,1);
    for i = 1:nt
        ftval(i) = cases(c).fTraw(Tvec(i));
    end
    ftval_norm = ftval / max(max(ftval), eps);

    % ODE simulation through the same shared right-hand side used by
    % sensitivity, elasticity and robustness analyses.
    p_case = p_model;
    p_case.T_opt = Topt_case;
    p_case.fT_norm = max(max(ftval), eps);  % identical observed-record normalization
    [E, N, A] = simulate_JL(Tvec, Hvec, Rvec, Wvec, dates, ht, p_case);

    % Validation against 2002 field data -> Figures 4-7 + Table 6 row
    [ValidationRow, kN] = validate_2002( ...
        thermoName, anno_cmp, dates, N, A, ...
        obsDatesA, obsAdultsTrap, obsDatesN, obsNymphsAllLeaf, ...
        fsAxes, fsTitle, fsLabel, fsLegend, saveFigures, figFolder);
    AllValidation = [AllValidation; ValidationRow]; %#ok<AGROW>

    % ENEI threshold calibration and rolling hindcast (Tables 7-8 data)
    [ThetaYear, RollingDiag, CurvesTable] = enei_calibration_LMH( ...
        thermoName, dates, E, N, Hvec, Rvec, ftval_norm, fH, fR, gammaE_max, ht, ...
        kN, years_eval, riskNames, riskVals, lead_days, ...
        thresholdStat, minPrevYears);
    AllThetaYear = [AllThetaYear; ThetaYear]; %#ok<AGROW>
    AllRolling   = [AllRolling;   RollingDiag]; %#ok<AGROW>
    AllCurves    = [AllCurves;    CurvesTable]; %#ok<AGROW>

    % Figure 9: calibration curves theta_r^(L) vs L
    plot_theta_vs_L(ThetaYear, thermoName, riskNames, lead_days, ...
        fsAxes, fsTitle, fsLabel, saveFigures, figFolder);

    % Figure 8: operational ENEI warning vs nymph risk crossings
    if thermoName == "Analytis"
        plot_warning_vs_risk_figure(CurvesTable, RollingDiag, thermoName, ...
            anno_cmp, selectedLeadForFigures, riskNames, riskVals, obs_tr_2002, ...
            fsAxes, fsTitle, fsLabel, fsLegend, saveFigures, figFolder);
    end
end

%% 12. PAPER TABLES: VALIDATION AND ENEI THRESHOLDS
% --- Table 6: validation metrics ---
T6 = table();
for c = 1:numel(cases)
    nm = string(cases(c).name);
    row = AllValidation(string(AllValidation.Thermal)==nm,:);
    T6 = [T6; table(nm, "Nymphs", row.kN(1), row.R2_Nymphs(1), row.p_Nymphs(1), ...
        row.AIC_Nymphs(1), row.RMSE_Nymphs(1), row.MAE_Nymphs(1), ...
        'VariableNames',{'Thermal','Stage','ScalingFactor','R2','p_value','AIC','RMSE','MAE'})]; %#ok<AGROW>
    T6 = [T6; table(nm, "Adults", row.kA(1), row.R2_Adults(1), row.p_Adults(1), ...
        row.AIC_Adults(1), row.RMSE_Adults(1), row.MAE_Adults(1), ...
        'VariableNames',{'Thermal','Stage','ScalingFactor','R2','p_value','AIC','RMSE','MAE'})]; %#ok<AGROW>
end
if saveTables
    writetable(T6, fullfile(tabFolder,'table06_validation_metrics.csv'));
end
fprintf('\n--- Table 6 (validation metrics) ---\n'); disp(T6)

% --- Table 7: rolling hindcast 2002 (low/medium/high x L=7/14/21) ---
T7 = AllRolling(AllRolling.Year==anno_cmp & ismember(AllRolling.TargetLeadDays,lead_days), ...
    {'Thermal','RiskLevel','TargetLeadDays','RealizedLeadDays','LeadErrorDays'});
T7 = sortrows(T7, {'Thermal','RiskLevel','TargetLeadDays'});
if saveTables
    writetable(T7, fullfile(tabFolder,'table07_rolling_hindcast_2002.csv'));
end
fprintf('\n--- Table 7 (rolling hindcast 2002) ---\n'); disp(T7)

% --- Table 8: median calibrated thresholds over 2000-2025 ---
leadCols = arrayfun(@(L) sprintf('L%d_days', L), lead_days, 'UniformOutput', false);
T8 = table();
for c = 1:numel(cases)
    nm = string(cases(c).name);
    for ir = 1:numel(riskNames)
        meds = nan(1,numel(lead_days));
        for il = 1:numel(lead_days)
            v = AllThetaYear.YearSpecificTheta( ...
                string(AllThetaYear.Thermal)==nm & ...
                AllThetaYear.RiskLevel==riskNames(ir) & ...
                AllThetaYear.TargetLeadDays==lead_days(il) & ...
                AllThetaYear.ValidTheta);
            meds(il) = median(v(isfinite(v)),'omitnan');
        end
        T8 = [T8; cell2table([{char(nm)} {char(riskNames(ir))} num2cell(meds)], ...
            'VariableNames',[{'Thermal','RiskLevel'} leadCols])]; %#ok<AGROW>
    end
end
if saveTables
    writetable(T8, fullfile(tabFolder,'table08_median_thresholds.csv'));
end
fprintf('\n--- Table 8 (median thresholds 2000-2025) ---\n'); disp(T8)

%% 13. MOVING-WINDOW OPERATIONAL THRESHOLDS AND LEAD-TIME BOXPLOT
% Analytis, W = 2 years, L = 21 days.
formulation_movwin = "Analytis";
TyA = AllThetaYear(string(AllThetaYear.Thermal)==formulation_movwin & ...
                   AllThetaYear.TargetLeadDays==selectedLeadForFigures & ...
                   AllThetaYear.ValidTheta, :);
CrA = AllCurves(string(AllCurves.Thermal)==formulation_movwin, :);

years_targ = 2002:2025;
nyr = numel(years_targ);
nrk = numel(riskNames);
theta_W2  = nan(nyr,nrk);
lead_W2   = nan(nyr,nrk);

for ii = 1:nyr
    yy = years_targ(ii);
    Cyy = CrA(CrA.Year==yy,:);
    if isempty(Cyy), continue; end
    datesY  = Cyy.Date; ENEI_y = Cyy.ENEI; Nsc_y = Cyy.Nscaled;
    for ir = 1:nrk
        prevYr = (yy-W_movwin):(yy-1);
        mPrev = ismember(TyA.Year, prevYr) & TyA.RiskLevel==riskNames(ir);
        prev = TyA.YearSpecificTheta(mPrev);
        prev = prev(isfinite(prev));
        if isempty(prev), continue; end
        thAgg = median(prev,'omitnan');
        theta_W2(ii,ir) = thAgg;
        idxR = find(Nsc_y >= riskVals(ir), 1, 'first');
        if isempty(idxR), continue; end
        tR   = datesY(idxR);
        idxW = find(ENEI_y >= thAgg, 1, 'first');
        if isempty(idxW), continue; end
        lead_W2(ii,ir) = days(tR - datesY(idxW));
    end
end

% Full-period reference (Table 8 medians) for the target L
ref_T8 = nan(1,nrk);
for ir = 1:nrk
    v = AllThetaYear.YearSpecificTheta( ...
        string(AllThetaYear.Thermal)==formulation_movwin & ...
        AllThetaYear.RiskLevel==riskNames(ir) & ...
        AllThetaYear.TargetLeadDays==selectedLeadForFigures & ...
        AllThetaYear.ValidTheta);
    ref_T8(ir) = median(v(isfinite(v)),'omitnan');
end

% --- Figure 10: year-by-year theta time series ---
% Colour map must have one row for each risk level.
% Low = blue, Medium = orange, High = purple.
if nrk == 3
    cols = [0.0000 0.4470 0.7410;   % Low
            0.8500 0.3250 0.0980;   % Medium
            0.4940 0.1840 0.5560];  % High
elseif nrk == 2
    cols = [0.8500 0.3250 0.0980;   % Medium
            0.4940 0.1840 0.5560];  % High
else
    cols = lines(nrk);
end
fig = figure('Color','w','Position',[100 100 1100 430]);
hold on
hLines = gobjects(nrk,1);
for ir = 1:nrk
    hLines(ir) = plot(years_targ, theta_W2(:,ir), 'o-', 'Color', cols(ir,:), ...
        'LineWidth', 1.4, 'MarkerFaceColor', cols(ir,:), 'MarkerSize', 5, ...
        'DisplayName', sprintf('%s (W=%d)', char(riskNames(ir)), W_movwin));
    m = isfinite(theta_W2(:,ir));
    if sum(m)>=3
        p = polyfit(years_targ(m), theta_W2(m,ir)', 1);
        plot(years_targ, polyval(p, years_targ), '--', 'Color', cols(ir,:), ...
            'LineWidth', 1.0, 'HandleVisibility','off');
    end
    yline(ref_T8(ir), ':', 'Color', cols(ir,:), 'LineWidth', 1.4, 'HandleVisibility','off');
    text(years_targ(end)+0.2, ref_T8(ir), ...
        sprintf('$\\widehat{\\theta}_{\\mathrm{%s}}^{(21)}$', ...
                char(riskNames(ir))), ...
        'Color', cols(ir,:), 'FontSize', 11, 'VerticalAlignment','middle', ...
        'Interpreter','latex');
end
set(gca,'YScale','log','XTick',2002:2:2025)
grid on
xlabel('Year','FontSize',fsLabel)
ylabel('$\widehat{\theta}^{(21)}_{r,y;\,W=2}$ (log scale)', ...
    'FontSize',fsLabel,'Interpreter','latex')
title('Year-by-year operational ENEI thresholds vs. full-period reference', ...
    'FontSize',fsTitle)
legend(hLines, 'Location','southeast','NumColumns',3,'FontSize',fsLegend)
set(gca,'FontSize',fsAxes)
hold off
if saveFigures
    disable_axes_toolbar(fig);
    exportgraphics(fig, fullfile(figFolder,'fig10_theta_movwin_W2.png'),'Resolution',300);
    exportgraphics(fig, fullfile(figFolder,'fig10_theta_movwin_W2.pdf'));
end

% --- Figure 11: boxplot of realised leads ---
fig = figure('Color','w','Position',[100 100 700 460]);
boxLabels = arrayfun(@char,riskNames,'UniformOutput',false);
boxplot(lead_W2, 'Labels', boxLabels, 'Widths', 0.55, ...
    'Symbol','o', 'OutlierSize',5);
hold on
% colour the boxes (Matlab returns Box patches in reverse order)
boxObjs = findobj(gca,'Tag','Box');
for k = 1:length(boxObjs)
    patch(get(boxObjs(k),'XData'), get(boxObjs(k),'YData'), ...
        cols(end-k+1,:), 'FaceAlpha',0.55, 'EdgeColor', cols(end-k+1,:));
end
% jittered points per year
rng(0);
for ir = 1:nrk
    v = lead_W2(:,ir); v = v(isfinite(v));
    xj = ir + (rand(numel(v),1)-0.5)*0.18;
    plot(xj, v, 'o', 'Color',[0.1 0.1 0.1], 'MarkerSize',3, 'MarkerFaceColor','none');
end
yline(selectedLeadForFigures, 'k:', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('target L = %d d', selectedLeadForFigures));
grid on
xlabel('Risk level','FontSize',fsLabel)
ylabel('Realised lead time (days)','FontSize',fsLabel)
title(sprintf('Distribution of realised lead times, 2002-2025 (W=%d, L=%d)', ...
    W_movwin, selectedLeadForFigures),'FontSize',fsTitle)
set(gca,'FontSize',fsAxes)
hold off
if saveFigures
    disable_axes_toolbar(fig);
    exportgraphics(fig, fullfile(figFolder,'fig11_leads_boxplot_W2.png'),'Resolution',300);
    exportgraphics(fig, fullfile(figFolder,'fig11_leads_boxplot_W2.pdf'));
end

% Distributional stats and yearly table export
fprintf('\n--- Realised lead distribution (W=%d, L=%d) ---\n', ...
    W_movwin, selectedLeadForFigures);
for ir = 1:nrk
    v = lead_W2(:,ir); v = v(isfinite(v));
    fprintf('  %-7s : median=%2d  IQR=[%2d,%2d]  range=[%2d,%2d]  mean=%4.1f +/- %3.1f\n', ...
        char(riskNames(ir)), round(median(v)), round(prctile(v,25)), round(prctile(v,75)), ...
        min(v), max(v), mean(v), std(v));
end

T_MW = table();
for ii=1:nyr
    for ir=1:nrk
        T_MW = [T_MW; table(years_targ(ii), string(riskNames(ir)), ...
            theta_W2(ii,ir), lead_W2(ii,ir), ...
            'VariableNames',{'Year','RiskLevel','ThetaW2','RealizedLeadDays'})]; %#ok<AGROW>
    end
end
if saveTables
    writetable(T_MW, fullfile(tabFolder,'theta_movwin_W2_yearly.csv'));
end

%% 14. TREND STATISTICS FOR MOVING-WINDOW THRESHOLDS
% OLS regression of theta_W2 on the target year, per risk level, with
% non-parametric Spearman rank correlation as a robustness check.
T9 = table();
for ir = 1:nrk
    y = theta_W2(:,ir);
    x = years_targ(:);
    mok = isfinite(y);
    if sum(mok) >= 3
        mdl   = fitlm(x(mok), y(mok));
        slope = mdl.Coefficients.Estimate(2);
        R2    = mdl.Rsquared.Ordinary;
        p_ols = mdl.Coefficients.pValue(2);
        [rho_sp, p_sp] = corr(x(mok), y(mok), 'Type', 'Spearman');
    else
        slope = NaN; R2 = NaN; p_ols = NaN; rho_sp = NaN; p_sp = NaN;
    end
    T9 = [T9; table(riskNames(ir), slope, R2, p_ols, rho_sp, p_sp, ...
        'VariableNames', {'RiskLevel','Slope_per_year','R2','p_OLS', ...
                          'Spearman_rho','p_Spearman'})]; %#ok<AGROW>
end
if saveTables
    writetable(T9, fullfile(tabFolder,'table09_theta_trend_W2.csv'));
end
fprintf('\n--- Table 9 (trend statistics of theta_W2, L=%d) ---\n', ...
    selectedLeadForFigures);
disp(T9)

%% 15. DEGREE-DAY BENCHMARK VS ENEI USING OBSERVED 2002 RISK DATES
% DD and Analytis-ENEI are evaluated against the same observed t_r.
% Degree-day and ENEI rolling hindcast for the target year 2002.
% The realised lead times are evaluated against interpolated observed
% nymph-risk crossing dates, not against model-generated Nscaled crossings.
% Calibration still uses the same model-based protocol for the preceding
% seasons, where observed nymph trajectories are not available.

formulation_dd_ref = "Analytis";

DD_year = nan(numel(years_eval), nrk, numel(lead_days));  % eta_{r,y}^{(L)}
DD_target_dates = NaT(0,0);
DD_target_traj  = [];

for iy = 1:numel(years_eval)
    yy = years_eval(iy);
    Cyy = AllCurves(string(AllCurves.Thermal)==formulation_dd_ref & ...
                    AllCurves.Year==yy, :);
    if isempty(Cyy), continue; end

    datesY = Cyy.Date;
    Ty     = Tvec(year(dates) == yy);
    Nsc_y  = Cyy.Nscaled;   % used for calibration years

    % biofix: first day of the year with T >= Tmin_spost
    bf = find(Ty >= Tmin_spost, 1, 'first');
    if isempty(bf), continue; end

    % daily DD increment with lower threshold Tmin and upper cut-off Tmax
    incDD = zeros(numel(Ty), 1);
    incDD(bf:end) = max(0, min(Ty(bf:end), Tmax) - Tmin);
    DDtr = cumsum(incDD);

    if yy == anno_cmp
        DD_target_dates = datesY;
        DD_target_traj  = DDtr;
    end

    % year-specific DD threshold values for each (risk, L)
    for ir = 1:nrk
        idxR = find(Nsc_y >= riskVals(ir), 1, 'first');
        if isempty(idxR), continue; end
        for il = 1:numel(lead_days)
            L = lead_days(il);
            iC = idxR - L;
            if iC >= 1
                DD_year(iy, ir, il) = DDtr(iC);
            end
        end
    end
end

% ENEI trajectory for the 2002 target season
C2002A = AllCurves(string(AllCurves.Thermal)=="Analytis" & ...
                   AllCurves.Year==anno_cmp, :);

T10 = table();
prev_iyy = find(ismember(years_eval, (anno_cmp - W_movwin):(anno_cmp - 1)));
prevYears = (anno_cmp - W_movwin):(anno_cmp - 1);

for ir = 1:nrk
    tR_obs = obs_tr_2002(ir);  % observed interpolated crossing date
    if isnat(tR_obs), continue; end

    for il = 1:numel(lead_days)
        L = lead_days(il);

        % --- DD branch ---
        prev_eta = squeeze(DD_year(prev_iyy, ir, il));
        prev_eta = prev_eta(isfinite(prev_eta));
        if isempty(prev_eta)
            eta_hat = NaN; twarn_dd = NaT; lead_dd = NaN; err_dd = NaN;
        else
            eta_hat = median(prev_eta, 'omitnan');
            idxWarn = find(DD_target_traj >= eta_hat, 1, 'first');
            if isempty(idxWarn)
                twarn_dd = NaT; lead_dd = NaN; err_dd = NaN;
            else
                twarn_dd = DD_target_dates(idxWarn);
                lead_dd  = days(tR_obs - twarn_dd);
                err_dd   = lead_dd - L;
            end
        end

        % --- ENEI Analytis branch ---
        maskPrevA = ismember(AllThetaYear.Year, prevYears) & ...
                    string(AllThetaYear.Thermal)=="Analytis" & ...
                    AllThetaYear.RiskLevel==riskNames(ir) & ...
                    AllThetaYear.TargetLeadDays==L & ...
                    AllThetaYear.ValidTheta;
        prevThetaA = AllThetaYear.YearSpecificTheta(maskPrevA);
        prevThetaA = prevThetaA(isfinite(prevThetaA));

        if isempty(prevThetaA) || isempty(C2002A)
            theta_hat_A = NaN; twarn_A = NaT; lead_A = NaN; err_A = NaN;
        else
            theta_hat_A = median(prevThetaA, 'omitnan');
            idxW = find(C2002A.ENEI >= theta_hat_A, 1, 'first');
            if isempty(idxW)
                twarn_A = NaT; lead_A = NaN; err_A = NaN;
            else
                twarn_A = C2002A.Date(idxW);
                lead_A  = days(tR_obs - twarn_A);
                err_A   = lead_A - L;
            end
        end

        T10 = [T10; table(riskNames(ir), L, string(tR_obs), ...
            eta_hat, string(twarn_dd), lead_dd, err_dd, ...
            theta_hat_A, string(twarn_A), lead_A, err_A, ...
            'VariableNames', {'RiskLevel','L_days','ObservedRiskDate', ...
                              'DD_eta_hat','DD_t_warn','DD_RealizedLead','DD_LeadError', ...
                              'ENEI_Analytis_theta','ENEI_Analytis_t_warn', ...
                              'ENEI_Analytis_Lead','ENEI_Analytis_LeadError'})]; %#ok<AGROW>
    end
end

T10 = sortrows(T10, {'RiskLevel','L_days'});
if saveTables
    writetable(T10, fullfile(tabFolder,'table10_dd_enei_vs_observed_tr.csv'));
    % Companion export using the previous Table 7 name, since the paper
    % comparison table is now based on the observed 2002 crossings.
    writetable(T10, fullfile(tabFolder,'table07_dd_enei_vs_observed_tr.csv'));
end
fprintf('\n--- Table 10 (DD vs ENEI, observed t_r, %d) ---\n', anno_cmp);
disp(T10)

fprintf('\nDone.\n');
fprintf('Figures saved in: %s\n', figFolder);
fprintf('Tables  saved in: %s\n', tabFolder);

%% 16. LOCAL FUNCTIONS
% The functions below keep the main script readable:
%   simulate_JL                         ODE time stepping
%   validate_2002                       scaling and validation metrics
%   plot_validation_regressions_2002    regression plots and statistics
%   enei_calibration_LMH                ENEI thresholds and rolling hindcast
%   plot_theta_vs_L                     calibration curves
%   plot_warning_vs_risk_figure         operational warning figure
%   fT_Analytis_raw                     thermal response function
%   disable_axes_toolbar                export helper

% ------------------------------------------------------------------------
% Simulate the stage-structured ODE model for one complete multi-year record.
% E = eggs, N = nymphs, A = adults, all expressed per leaf before scaling.
function [E,N,A] = simulate_JL(Tvec, Hvec, Rvec, Wvec, dates, ht, p)
% SIMULATE_JL  Multi-year integration wrapper around the shared enei_rhs.
% Annual reset and immigration logic are kept here; every biological rate and
% every environmental term is evaluated in enei_rhs.m.

    nt = length(Tvec);
    E = zeros(nt,1);
    N = zeros(nt,1);
    A = zeros(nt,1);
    adults_introduced_this_year = false;

    for j = 1:nt-1
        doy = day(dates(j),'dayofyear');

        if doy == 1
            E(j)=0; N(j)=0; A(j)=0;
            adults_introduced_this_year = false;
        end

        if (Tvec(j) > p.T_move) && ~adults_introduced_this_year
            A(j) = A(j) + p.A_intro;
            adults_introduced_this_year = true;
        end

        [dstate, ~] = enei_rhs([E(j); N(j); A(j)], ...
            Tvec(j), Hvec(j), Rvec(j), Wvec(j), p);

        E(j+1) = max(0, E(j) + ht*dstate(1));
        N(j+1) = max(0, N(j) + ht*dstate(2));
        A(j+1) = max(0, A(j) + ht*dstate(3));
    end
end

% ------------------------------------------------------------------------
% Compare the 2002 model trajectory with field observations.
% Peak-based kA and kN scale model outputs to observation units.
function [ValidationRow, kN] = validate_2002( ...
    thermoName, anno_cmp, dates, N, A, obsDatesA, obsAdultsTrap, obsDatesN, obsNymphsAllLeaf, ...
    fsAxes, fsTitle, fsLabel, fsLegend, saveFigures, figFolder)

    thermoName = string(thermoName);
    thermoChar = char(thermoName);

    maskY = year(dates) == anno_cmp;
    datesY = dates(maskY); datesY = datesY(:);
    NY = N(maskY); NY = NY(:);
    AY = A(maskY); AY = AY(:);

    A_mod_obs_raw = interp1(datenum(datesY), AY, datenum(obsDatesA), 'linear', 'extrap');
    N_mod_obs_raw = interp1(datenum(datesY), NY, datenum(obsDatesN), 'linear', 'extrap');

    [obsPeakA, iObsPeakA] = max(obsAdultsTrap);
    [modPeakA, iModPeakA] = max(AY);
    [obsPeakN, iObsPeakN] = max(obsNymphsAllLeaf);
    [modPeakN, iModPeakN] = max(NY);

    kA = modPeakA / obsPeakA;
    kN = modPeakN / obsPeakN;

    AY_scaled = AY / kA;
    NY_scaled = NY / kN;
    A_mod_obs_scaled = A_mod_obs_raw / kA;
    N_mod_obs_scaled = N_mod_obs_raw / kN;

    rmseA = sqrt(mean((A_mod_obs_scaled - obsAdultsTrap).^2, 'omitnan'));
    maeA  = mean(abs(A_mod_obs_scaled - obsAdultsTrap), 'omitnan');
    rmseN = sqrt(mean((N_mod_obs_scaled - obsNymphsAllLeaf).^2, 'omitnan'));
    maeN  = mean(abs(N_mod_obs_scaled - obsNymphsAllLeaf), 'omitnan');

    [peakA_scaled, iPeakA_scaled] = max(AY_scaled);
    [peakN_scaled, iPeakN_scaled] = max(NY_scaled);

    modelAdultPeakDate = datesY(iPeakA_scaled);
    modelNymphPeakDate = datesY(iPeakN_scaled);
    obsAdultPeakDate = obsDatesA(iObsPeakA);
    obsNymphPeakDate = obsDatesN(iObsPeakN);

    fprintf('\n[%s] Adults scaling: kA = %.4f | RMSE = %.3f | MAE = %.3f\n', thermoChar, kA, rmseA, maeA);
    fprintf('[%s] Nymphs scaling: kN = %.4f | RMSE = %.3f | MAE = %.3f\n', thermoChar, kN, rmseN, maeN);

    % Time-series validation figure
    fig = figure('Color','w','Position',[100 100 1100 430]);
    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile
    plot(datesY, NY_scaled, 'LineWidth', 2, 'DisplayName','Model scaled'); hold on
    scatter(obsDatesN, obsNymphsAllLeaf, 45, 'filled', 'DisplayName','Observed'); hold off
    grid on
    xlabel('Date','FontSize',fsLabel)
    ylabel('Nymphs per leaf','FontSize',fsLabel)
    title(sprintf('Nymphs, %d', anno_cmp), 'FontSize',fsTitle)
    legend('Location','best','FontSize',fsLegend)
    set(gca,'FontSize',fsAxes)

    nexttile
    plot(datesY, AY_scaled, 'LineWidth', 2, 'DisplayName','Model scaled'); hold on
    scatter(obsDatesA, obsAdultsTrap, 45, 'filled', 'DisplayName','Observed'); hold off
    grid on
    xlabel('Date','FontSize',fsLabel)
    ylabel('Adults per trap','FontSize',fsLabel)
    title(sprintf('Adults, %d', anno_cmp), 'FontSize',fsTitle)
    legend('Location','best','FontSize',fsLegend)
    set(gca,'FontSize',fsAxes)

    title(tl, sprintf('Feudo Arancio %d validation', anno_cmp), 'FontSize', fsTitle)

    if saveFigures
        disable_axes_toolbar(fig);
        exportgraphics(fig, fullfile(figFolder, sprintf('fig_validation_timeseries_%s_%d.png', thermoChar, anno_cmp)), 'Resolution', 300);
        exportgraphics(fig, fullfile(figFolder, sprintf('fig_validation_timeseries_%s_%d.pdf', thermoChar, anno_cmp)));
    end

    % Regression validation figure
    RegStats = plot_validation_regressions_2002( ...
        thermoName, anno_cmp, obsNymphsAllLeaf, N_mod_obs_scaled, ...
        obsAdultsTrap, A_mod_obs_scaled, fsAxes, fsTitle, fsLabel, fsLegend, ...
        saveFigures, figFolder);

    ValidationRow = table( ...
        thermoName, kA, kN, rmseA, maeA, rmseN, maeN, ...
        obsPeakA, modPeakA, peakA_scaled, string(obsAdultPeakDate), string(datesY(iModPeakA)), string(modelAdultPeakDate), ...
        obsPeakN, modPeakN, peakN_scaled, string(obsNymphPeakDate), string(datesY(iModPeakN)), string(modelNymphPeakDate), ...
        RegStats.R2_Nymphs, RegStats.p_Nymphs, RegStats.AIC_Nymphs, ...
        RegStats.R2_Adults, RegStats.p_Adults, RegStats.AIC_Adults, ...
        'VariableNames', { ...
        'Thermal','kA','kN','RMSE_Adults','MAE_Adults','RMSE_Nymphs','MAE_Nymphs', ...
        'ObservedAdultPeak','RawModelAdultPeak','ScaledModelAdultPeak', ...
        'ObservedAdultPeakDate','RawModelAdultPeakDate','ScaledModelAdultPeakDate', ...
        'ObservedNymphPeak','RawModelNymphPeak','ScaledModelNymphPeak', ...
        'ObservedNymphPeakDate','RawModelNymphPeakDate','ScaledModelNymphPeakDate', ...
        'R2_Nymphs','p_Nymphs','AIC_Nymphs','R2_Adults','p_Adults','AIC_Adults'});
end

% ------------------------------------------------------------------------
% Build observed-vs-modelled regression plots and return R2, p-values and AIC.
function RegStats = plot_validation_regressions_2002( ...
    thermoName, anno_cmp, obsNymph, modNymph, obsAdult, modAdult, ...
    fsAxes, fsTitle, fsLabel, fsLegend, saveFigures, figFolder)

    thermoName = string(thermoName);
    thermoChar = char(thermoName);

    maskN = isfinite(obsNymph) & isfinite(modNymph);
    xN = obsNymph(maskN); yN = modNymph(maskN);

    maskA = isfinite(obsAdult) & isfinite(modAdult);
    xA = obsAdult(maskA); yA = modAdult(maskA);

    mdlN = fitlm(xN, yN);
    mdlA = fitlm(xA, yA);

    xgN = linspace(min(xN), max(xN), 150)';
    xgA = linspace(min(xA), max(xA), 150)';
    [yfitN, yciN] = predict(mdlN, xgN, 'Alpha', 0.05);
    [yfitA, yciA] = predict(mdlA, xgA, 'Alpha', 0.05);

    RegStats = struct();
    RegStats.R2_Nymphs  = mdlN.Rsquared.Ordinary;
    RegStats.p_Nymphs   = mdlN.Coefficients.pValue(2);
    RegStats.AIC_Nymphs = mdlN.ModelCriterion.AIC;
    RegStats.R2_Adults  = mdlA.Rsquared.Ordinary;
    RegStats.p_Adults   = mdlA.Coefficients.pValue(2);
    RegStats.AIC_Adults = mdlA.ModelCriterion.AIC;

    fig = figure('Color','w','Position',[100 100 1200 450]);
    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile
    hData = scatter(xN, yN, 45, '^', 'filled', 'DisplayName','Data'); hold on
    hFit  = plot(xgN, yfitN, 'k-', 'LineWidth', 1.6, 'DisplayName','Fit');
    hCI   = plot(xgN, yciN(:,1), 'r--', 'LineWidth', 1.2, 'DisplayName','95% CI');
    plot(xgN, yciN(:,2), 'r--', 'LineWidth', 1.2, 'HandleVisibility','off'); hold off
    grid on
    xlabel('Observed data (Nymph)', 'FontSize', fsLabel)
    ylabel('Modelled data (Nymph)', 'FontSize', fsLabel)
    title(sprintf('Feudo Arancio %d: nymphs linear regression', anno_cmp), 'FontSize', fsTitle)
    legend([hData hFit hCI], 'Location','southeast', 'FontSize',fsLegend)
    set(gca,'FontSize',fsAxes)
    txtN = sprintf('R^2 = %.2f\np = %.3f\nAIC = %.1f', RegStats.R2_Nymphs, RegStats.p_Nymphs, RegStats.AIC_Nymphs);
    text(0.03, 0.95, txtN, 'Units','normalized', 'VerticalAlignment','top', ...
        'BackgroundColor','w', 'EdgeColor',[0.7 0.7 0.7], 'FontSize',fsAxes);

    nexttile
    hData = scatter(xA, yA, 45, 'o', 'filled', 'DisplayName','Data'); hold on
    hFit  = plot(xgA, yfitA, 'k-', 'LineWidth', 1.6, 'DisplayName','Fit');
    hCI   = plot(xgA, yciA(:,1), 'r--', 'LineWidth', 1.2, 'DisplayName','95% CI');
    plot(xgA, yciA(:,2), 'r--', 'LineWidth', 1.2, 'HandleVisibility','off'); hold off
    grid on
    xlabel('Observed data (Adult)', 'FontSize', fsLabel)
    ylabel('Modelled data (Adult)', 'FontSize', fsLabel)
    title(sprintf('Feudo Arancio %d: adults linear regression', anno_cmp), 'FontSize', fsTitle)
    legend([hData hFit hCI], 'Location','southeast', 'FontSize',fsLegend)
    set(gca,'FontSize',fsAxes)
    txtA = sprintf('R^2 = %.2f\np = %.3f\nAIC = %.1f', RegStats.R2_Adults, RegStats.p_Adults, RegStats.AIC_Adults);
    text(0.03, 0.95, txtA, 'Units','normalized', 'VerticalAlignment','top', ...
        'BackgroundColor','w', 'EdgeColor',[0.7 0.7 0.7], 'FontSize',fsAxes);

    title(tl, 'Linear regression analysis', 'FontSize', fsTitle)

    if saveFigures
        disable_axes_toolbar(fig);
        exportgraphics(fig, fullfile(figFolder, sprintf('fig_regression_2002_%s.png', thermoChar)), 'Resolution', 300);
        exportgraphics(fig, fullfile(figFolder, sprintf('fig_regression_2002_%s.pdf', thermoChar)));
    end
end

% ------------------------------------------------------------------------
% Compute yearly ENEI curves, year-specific thresholds theta_r^(L),
% and rolling-hindcast warning diagnostics for Low/Medium/High risk levels.
function [ThetaYear, RollingDiag, CurvesTable] = enei_calibration_LMH( ...
    thermoName, dates, E, N, Hvec, Rvec, ftval_norm, fH, fR, gammaE_max, ht, ...
    kN_scale, years_eval, riskNames, riskVals, lead_days, thresholdStat, minPrevYears)

    thermoName = string(thermoName);
    thresholdStat = lower(string(thresholdStat));
    ThetaYear = table();
    RollingDiag = table();
    CurvesTable = table();

    Y = struct('Year',{},'dates',{},'ENEI',{},'Nscaled',{},'idxRisk',{},'tRisk',{});

    for yy = years_eval(:)'
        maskY = year(dates)==yy;
        if ~any(maskY), continue; end

        idxY = find(maskY);
        datesY = dates(maskY); datesY = datesY(:);
        EY = E(maskY); EY = EY(:);
        NYraw = N(maskY); NYraw = NYraw(:);
        HY = Hvec(idxY); HY = HY(:);
        RY = Rvec(idxY); RY = RY(:);
        ftY = ftval_norm(idxY); ftY = ftY(:);

        gammaE_Y = gammaE_max .* ftY .* fH(HY) .* fR(RY);
        ENEI = cumsum(gammaE_Y .* EY) * ht;
        Nscaled = NYraw ./ kN_scale;

        idxRisk = nan(1, numel(riskVals));
        tRisk = NaT(1, numel(riskVals));
        for ir = 1:numel(riskVals)
            idx = find(Nscaled >= riskVals(ir), 1, 'first');
            if ~isempty(idx)
                idxRisk(ir) = idx;
                tRisk(ir) = datesY(idx);
            end
        end

        CurvesTable = [CurvesTable; table( ...
            repmat(thermoName,numel(datesY),1), repmat(yy,numel(datesY),1), ...
            datesY, day(datesY,'dayofyear'), ENEI, Nscaled, ...
            'VariableNames', {'Thermal','Year','Date','DayOfYear','ENEI','Nscaled'})]; %#ok<AGROW>

        Y(end+1).Year = yy; %#ok<AGROW>
        Y(end).dates = datesY;
        Y(end).ENEI = ENEI;
        Y(end).Nscaled = Nscaled;
        Y(end).idxRisk = idxRisk;
        Y(end).tRisk = tRisk;

        for ir = 1:numel(riskVals)
            for il = 1:numel(lead_days)
                L = lead_days(il);
                lagSteps = round(L / ht);
                theta = NaN; calDate = NaT; valid = false;
                idxr = idxRisk(ir);
                if isfinite(idxr)
                    idxCal = idxr - lagSteps;
                    if idxCal >= 1
                        theta = ENEI(idxCal);
                        calDate = datesY(idxCal);
                        valid = isfinite(theta);
                    end
                end
                ThetaYear = [ThetaYear; table( ...
                    thermoName, yy, riskNames(ir), riskVals(ir), L, theta, calDate, tRisk(ir), valid, ...
                    'VariableNames', {'Thermal','Year','RiskLevel','RiskThreshold','TargetLeadDays', ...
                                      'YearSpecificTheta','CalibrationDate','RiskDate','ValidTheta'})]; %#ok<AGROW>
            end
        end
    end

    for iy = 1:numel(Y)
        yy = Y(iy).Year;
        datesY = Y(iy).dates;
        ENEI = Y(iy).ENEI;

        for ir = 1:numel(riskVals)
            tRisk = Y(iy).tRisk(ir);
            idxr = Y(iy).idxRisk(ir);

            for il = 1:numel(lead_days)
                L = lead_days(il);
                maskPrev = ThetaYear.Year < yy & ...
                           ThetaYear.RiskLevel == riskNames(ir) & ...
                           ThetaYear.TargetLeadDays == L & ...
                           ThetaYear.ValidTheta;
                prevVals = ThetaYear.YearSpecificTheta(maskPrev);
                prevVals = prevVals(isfinite(prevVals));
                nPrev = numel(prevVals);

                thetaAgg = NaN;
                if nPrev >= minPrevYears
                    switch thresholdStat
                        case "median"
                            thetaAgg = median(prevVals, 'omitnan');
                        case "q25"
                            thetaAgg = prctile(prevVals, 25);
                        case "q10"
                            thetaAgg = prctile(prevVals, 10);
                        otherwise
                            thetaAgg = median(prevVals, 'omitnan');
                    end
                end

                warningDate = NaT;
                thresholdReached = false;
                realizedLead = NaN;
                leadError = NaN;
                if isfinite(thetaAgg)
                    idxWarn = find(ENEI >= thetaAgg, 1, 'first');
                    if ~isempty(idxWarn)
                        warningDate = datesY(idxWarn);
                        thresholdReached = true;
                    end
                end
                if thresholdReached && isfinite(idxr)
                    realizedLead = days(tRisk - warningDate);
                    leadError = realizedLead - L;
                end

                RollingDiag = [RollingDiag; table( ...
                    thermoName, yy, riskNames(ir), riskVals(ir), L, thetaAgg, warningDate, tRisk, ...
                    realizedLead, leadError, nPrev, thresholdReached, ...
                    'VariableNames', {'Thermal','Year','RiskLevel','RiskThreshold','TargetLeadDays', ...
                                      'AggregatedTheta','WarningDate','RiskDate','RealizedLeadDays', ...
                                      'LeadErrorDays','NPreviousYears','ThresholdReached'})]; %#ok<AGROW>
            end
        end
    end
end

% ------------------------------------------------------------------------
% Plot median ENEI calibration curves theta_r^(L) versus warning horizon L.
function plot_theta_vs_L(ThetaYear, thermoName, riskNames, lead_days, ...
    fsAxes, fsTitle, fsLabel, saveFigures, figFolder)

    thermoChar = char(string(thermoName));
    D = ThetaYear(string(ThetaYear.Thermal)==string(thermoName) & ThetaYear.ValidTheta,:);
    if isempty(D), return; end

    fig = figure('Color','w','Position',[100 100 1200 320]);
    tl = tiledlayout(1,numel(riskNames),'TileSpacing','compact','Padding','compact');

    for ir = 1:numel(riskNames)
        nexttile
        medTheta = nan(size(lead_days));
        q25Theta = nan(size(lead_days));
        q75Theta = nan(size(lead_days));
        for il = 1:numel(lead_days)
            vals = D.YearSpecificTheta(D.RiskLevel==riskNames(ir) & D.TargetLeadDays==lead_days(il));
            vals = vals(isfinite(vals));
            if ~isempty(vals)
                medTheta(il)=median(vals,'omitnan');
                q25Theta(il)=prctile(vals,25);
                q75Theta(il)=prctile(vals,75);
            end
        end
        errorbar(lead_days, medTheta, medTheta-q25Theta, q75Theta-medTheta, '-o','LineWidth',1.8,'MarkerSize',7);
        finiteMed = medTheta(isfinite(medTheta));
        if ~isempty(finiteMed) && all(finiteMed > 0), set(gca,'YScale','log'); end
        grid on
        xticks(lead_days)
        xlabel('Lead time L [days]','FontSize',fsLabel)
        ylabel('\theta_r^{(L)}','FontSize',fsLabel)
        title(sprintf('%s risk', char(riskNames(ir))),'FontSize',fsTitle)
        set(gca,'FontSize',fsAxes)
    end
    title(tl, 'Calibration curves','FontSize',fsTitle)

    if saveFigures
        disable_axes_toolbar(fig);
        exportgraphics(fig, fullfile(figFolder, sprintf('fig01_theta_vs_L_%s.png', thermoChar)), 'Resolution', 300);
        exportgraphics(fig, fullfile(figFolder, sprintf('fig01_theta_vs_L_%s.pdf', thermoChar)));
    end
end

% ------------------------------------------------------------------------
% Plot the operational warning figure for the target year.
% The bottom panel uses observed risk-crossing dates when provided.
function plot_warning_vs_risk_figure(CurvesTable, RollingDiag, thermoName, ...
    anno_cmp, selectedLead, riskNames, riskVals, obs_tr_2002, ...
    fsAxes, fsTitle, fsLabel, fsLegend, saveFigures, figFolder)

    thermoName = string(thermoName);
    thermoChar = char(thermoName);
    riskNames = string(riskNames);
    nrk = numel(riskNames);

    % ===================== Load target-year data =========================
    C = CurvesTable(string(CurvesTable.Thermal)==thermoName & ...
                    CurvesTable.Year==anno_cmp, :);

    R = RollingDiag(string(RollingDiag.Thermal)==thermoName & ...
                    RollingDiag.Year==anno_cmp & ...
                    RollingDiag.TargetLeadDays==selectedLead, :);

    if isempty(C) || isempty(R)
        warning('No data available for %s, year %d, L=%d.', ...
            thermoChar, anno_cmp, selectedLead);
        return
    end

    % Keep only the risk levels requested by riskNames and sort them.
    keep = ismember(string(R.RiskLevel), riskNames);
    R = R(keep,:);
    [~, ord] = ismember(string(R.RiskLevel), riskNames);
    [~, idxSort] = sort(ord);
    R = R(idxSort,:);

    if height(R) ~= nrk
        warning('Missing one or more requested risk levels for %s, year %d.', ...
            thermoChar, anno_cmp);
        return
    end

    thetaWarn = R.AggregatedTheta(:);
    tWarn     = R.WarningDate(:);

    if nargin < 8 || isempty(obs_tr_2002)
        tRisk = R.RiskDate(:);
        riskDateLabel = 'model-based';
    else
        tRisk = obs_tr_2002(:);
        riskDateLabel = 'observed';
    end

    if numel(tRisk) ~= nrk
        warning('obs_tr_2002 has incompatible length for %s, year %d.', ...
            thermoChar, anno_cmp);
        return
    end

    realisedLead = days(tRisk - tWarn);

    if any(~isfinite(thetaWarn)) || any(isnat(tWarn)) || any(isnat(tRisk))
        warning('Incomplete warning/risk dates or thresholds for %s, year %d, L=%d.', ...
            thermoChar, anno_cmp, selectedLead);
        return
    end

    % ===================== Plot window ===================================
    validDates = [tWarn(:); tRisk(:)];
    validDates = validDates(~isnat(validDates));

    if isempty(validDates)
        x0 = min(C.Date);
        x1 = max(C.Date);
    else
        x0 = min(validDates) - days(10);
        x1 = max(validDates) + days(10);
    end

    maskWin = C.Date >= x0 & C.Date <= x1;
    Cw = C(maskWin,:);

    % ===================== Colors ========================================
    if nrk == 2
        colsRisk = [0.8500 0.3250 0.0980;   % medium: orange
                    0.4940 0.1840 0.5560];  % high: purple
    elseif nrk == 3
        colsRisk = [0.0000 0.4470 0.7410;   % low: blue
                    0.8500 0.3250 0.0980;   % medium: orange
                    0.4940 0.1840 0.5560];  % high: purple
    else
        colsRisk = lines(nrk);
    end
    colNoWarning = [0.65 0.65 0.65];
    colRiskCurve = [0.2 0.2 0.2];

    % ===================== Figure ========================================
    fig = figure('Color','w','Position',[100 100 1250 950]);
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    % =====================================================================
    % PANEL 1: ENEI(t) and warning thresholds
    % =====================================================================
    nexttile
    hold on
    grid on

    plot(Cw.Date, Cw.ENEI, 'k-', 'LineWidth', 2.0, 'DisplayName','ENEI(t)');

    for ir = 1:nrk
        yline(thetaWarn(ir), '--', 'Color', colsRisk(ir,:), 'LineWidth', 1.4, ...
            'DisplayName', sprintf('$\\widehat{\\theta}_{\\mathrm{%s}}^{(L)}$', ...
            lower(char(riskNames(ir)))));
        xline(tWarn(ir), '--', 'Color', colsRisk(ir,:), 'LineWidth', 1.3, ...
            'HandleVisibility','off');
        plot(tWarn(ir), thetaWarn(ir), 'o', 'Color', colsRisk(ir,:), ...
            'MarkerFaceColor','w', 'MarkerSize', 7, 'HandleVisibility','off');
        text(tWarn(ir), thetaWarn(ir), ...
            sprintf('  %s warning\n  %s', lower(char(riskNames(ir))), ...
            datestr(tWarn(ir),'dd-mmm')), ...
            'Color', colsRisk(ir,:), 'VerticalAlignment','bottom');
    end

    xlim([x0 x1])
    xlabel('Date','FontSize',fsLabel)
    ylabel('ENEI (log scale)','FontSize',fsLabel)
    title(sprintf('ENEI warning thresholds, %d, L=%d days', ...
        anno_cmp, selectedLead), 'FontSize', fsTitle)
    legend('Location','northwest','FontSize',fsLegend,'Interpreter','latex')
    set(gca,'FontSize',fsAxes,'YScale','log')

    positiveVals = [thetaWarn(:); Cw.ENEI(Cw.ENEI > 0)];
    if isempty(positiveVals)
        ylim([1 10])
    else
        yl = ylim;
        ylim([max(eps, min(positiveVals)*0.3), max(yl(2), max(positiveVals)*1.15)])
    end

    hold off

    % =====================================================================
    % PANEL 2: ENEI warning calendar
    % =====================================================================
    nexttile
    hold on
    grid on

    yBar = 1.0;
    lwBar = 8;

    % segments based on warning thresholds crossed by ENEI
    % no warning | low warning | medium warning | high warning (for the current settings)
    plot([x0 tWarn(1)], [yBar yBar], '-', 'Color', colNoWarning, ...
        'LineWidth', lwBar);

    for ir = 1:nrk
        if ir < nrk
            segEnd = tWarn(ir+1);
        else
            segEnd = x1;
        end
        plot([tWarn(ir) segEnd], [yBar yBar], '-', 'Color', colsRisk(ir,:), ...
            'LineWidth', lwBar);
        xline(tWarn(ir), '--', 'Color', colsRisk(ir,:), 'LineWidth', 1.2, ...
            'HandleVisibility','off');
        text(tWarn(ir), 1.13, ...
            sprintf('%s warning\n%s', lower(char(riskNames(ir))), ...
            datestr(tWarn(ir),'dd-mmm')), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'Color', colsRisk(ir,:));
    end

    text(x0 + (tWarn(1)-x0)/2, 0.82, 'no warning', ...
        'HorizontalAlignment','center', 'Color', [0.4 0.4 0.4]);

    for ir = 1:nrk
        if ir < nrk
            segEnd = tWarn(ir+1);
        else
            segEnd = x1;
        end
        text(tWarn(ir) + (segEnd-tWarn(ir))/2, 0.82, ...
            sprintf('%s warning state', lower(char(riskNames(ir)))), ...
            'HorizontalAlignment','center', 'Color', colsRisk(ir,:));
    end

    xlim([x0 x1])
    ylim([0.65 1.25])
    yticks([])
    xlabel('Date','FontSize',fsLabel)
    title('ENEI warning calendar', 'FontSize', fsTitle)
    set(gca,'FontSize',fsAxes)

    hold off

    % =====================================================================
    % PANEL 3: Nscaled(t) and observed risk crossings
    % =====================================================================
    nexttile
    hold on
    grid on

    hN = plot(Cw.Date, Cw.Nscaled, 'Color', colRiskCurve, 'LineWidth', 2.0, ...
        'DisplayName','N_{scaled}(t)');

    hRisk = gobjects(nrk,1);
    for ir = 1:nrk
        hRisk(ir) = yline(riskVals(ir), '--', 'Color', colsRisk(ir,:), ...
            'LineWidth', 1.3, ...
            'DisplayName', sprintf('N_{%s}', lower(char(riskNames(ir)))));
        xline(tRisk(ir), '-', 'Color', colsRisk(ir,:), 'LineWidth', 1.3, ...
            'HandleVisibility','off');
        plot(tRisk(ir), riskVals(ir), 'o', 'Color', colsRisk(ir,:), ...
            'MarkerFaceColor','w', 'MarkerSize', 7, 'HandleVisibility','off');
        text(tRisk(ir), riskVals(ir), ...
            sprintf('  %s risk\n  %s', lower(char(riskNames(ir))), ...
            datestr(tRisk(ir),'dd-mmm')), ...
            'Color', colsRisk(ir,:), 'VerticalAlignment','bottom');
    end

    % lead-time arrows (warning date -> observed risk-crossing date)
    yTop = max([Cw.Nscaled; riskVals(:)]) * 1.12;
    if isempty(yTop) || ~isfinite(yTop) || yTop <= 0
        yTop = max(riskVals) * 1.12;
    end

    for ir = 1:nrk
        yy = yTop * (1 - 0.07*(ir-1));
        plot([tWarn(ir) tRisk(ir)], [yy yy], '-', 'Color', colsRisk(ir,:), ...
            'LineWidth', 1.5, 'HandleVisibility','off');
        plot(tWarn(ir), yy, '<', 'Color', colsRisk(ir,:), ...
            'MarkerFaceColor', colsRisk(ir,:), 'HandleVisibility','off');
        plot(tRisk(ir), yy, '>', 'Color', colsRisk(ir,:), ...
            'MarkerFaceColor', colsRisk(ir,:), 'HandleVisibility','off');
        text(tWarn(ir) + (tRisk(ir)-tWarn(ir))/2, yy*1.02, ...
            sprintf('%.0f d', realisedLead(ir)), ...
            'HorizontalAlignment','center', 'Color', colsRisk(ir,:));
    end

    ylim([0, yTop*1.15])
    xlim([x0 x1])
    xlabel('Date','FontSize',fsLabel)
    ylabel('N_{scaled}(t) [nymphs/leaf]','FontSize',fsLabel)
    title(sprintf('Scaled nymph density and %s risk crossings', riskDateLabel), ...
        'FontSize', fsTitle)
    legend([hN; hRisk], 'Location','northwest', 'FontSize', fsLegend)
    set(gca,'FontSize',fsAxes)

    hold off

    % =====================================================================
    % Overall title
    % =====================================================================
    title(tl, sprintf('ENEI warning states versus observed nymph-based risk crossings (%d, L=%d)', ...
        anno_cmp, selectedLead), 'FontSize', fsTitle+1);

    % ===================== Save ==========================================
    if saveFigures
        disable_axes_toolbar(fig);
        exportgraphics(fig, fullfile(figFolder, ...
            sprintf('fig_warning_vs_observed_risk_%s_%d_L%d.png', ...
            thermoChar, anno_cmp, selectedLead)), 'Resolution', 300);

        exportgraphics(fig, fullfile(figFolder, ...
            sprintf('fig_warning_vs_observed_risk_%s_%d_L%d.pdf', ...
            thermoChar, anno_cmp, selectedLead)));
    end
end

% ------------------------------------------------------------------------

function f = fT_Analytis_raw(T, Tmin, Tmax)
    if T <= Tmin || T >= Tmax
        f = 0;
    else
        f = (T - Tmin) * sqrt(max(Tmax - T, 0));
    end
end

% ------------------------------------------------------------------------
function disable_axes_toolbar(fig)
    % Avoid exporting the interactive axes toolbar in saved figures.
    try
        ax = findall(fig, 'Type', 'Axes');
        for ia = 1:numel(ax)
            try
                ax(ia).Toolbar.Visible = 'off';
            catch
            end
        end
    catch
    end
end
