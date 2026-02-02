% ==============================================================
% Jacobiasca lybica Stage-Structured Model (per leaf)
% Two thermal cases:
%   1) Analytis:  f_T = (T-Tmin)*sqrt(Tmax-T)
%   2) Briere:    f_T = T*(T-Tmin)*sqrt(Tmax-T)
%
% - Feudo Arancio 2002 comparison using peak-based kA and kN
%   (scales the MODEL to the units of the observations)
% - Early warning 
% ==============================================================

clc; clear; close all

%% ===================== METEOROLOGICAL DATA ===============================
load('Arancio_dati_2000_2025.mat')   % Tvec, Hvec, Rvec, Wvec

%% ===================== BIOLOGICAL PARAMETERS ======================
% Reproduction and eggs
beta_max   = 5;        % [eggs d^-1 adult^-1]
muE_min    = 0.05;     % [d^-1]
alphaT_E   = 0.02;     % [d^-1 °C^-1]
alphaH_E   = 0.003;    % [d^-1 %^-1]
gammaE_max = 0.235;    % [d^-1] max tasso schiusa

% Nymphs
muN_min    = 0.042;    % [d^-1]
kW_mort    = 0.10;     % [d^-1]
gammaN_max = 0.12;     % [d^-1]
alphaT_N   = 0.006;    % [d^-1 °C^-1]
alphaH_N   = 0.0012;   % [d^-1 %^-1]

% Adulti
muA_min  = 0.03;       % [d^-1]
betaT_A  = 0.0145;     % [d^-1 °C^-1]
betaH_A  = 0.002;      % [d^-1 %^-1]

%% ===================== ENVIRONMENTAL FACTORS =======================
Tmin       = 11.0;     % [°C]
Tmax       = 35.0;     % [°C]
Tmin_spost = 22;       % [°C] adult introduction threshold (T_move)

Topt_HR = 60.0;        % [%] optimal relative humidity
Hopt    = 60.0;        % [%]
cH      = 0.08;        % logistic slope for H
H0      = 60.0;        % logistic midpoint (center) for H

cR = 0.25;             % logistic slope for R
R0 = 18.0;             % logistic midpoint (center) for R

KW_half    = 3.0;      % [m s^-1] half-saturation wind speed
alphaW_dev = 0.10;     % [-] dimensionless


%% ===================== NON-THERMAL FUNCTIONS ====================
fH = @(H) 1.0 ./ (1.0 + exp(-cH*(H - H0)));
fR = @(R) 1.0 ./ (1.0 + exp(-cR*(R - R0)));
fW = @(W) W ./ (W + KW_half);
gW = @(W) 1 ./ (1 + alphaW_dev*W);

%% ===================== TIME ===================================
nt = length(Tvec);
startDate = datetime(2000,1,1);
dates     = startDate + days(0:nt-1);
ht        = 1;

%% ===================== THERMAL CASES =============================
cases(1).name  = "Analytis";
cases(1).fTraw = @(T) fT_Analytis_raw(T, Tmin, Tmax);
cases(1).Topt  = (2*Tmax + Tmin)/3;   % n=1, m=0.5

cases(2).name  = "Briere";
cases(2).fTraw = @(T) fT_Briere_raw(T, Tmin, Tmax);
cases(2).Topt  = briere_Topt_classic(Tmin, Tmax);

%% ===================== SIMULATIONS (TWO CASES) ===================
Results = struct();

for c = 1:numel(cases)
    thermoName = cases(c).name;
    Topt_case  = cases(c).Topt;

    % ---- % Daily f_T computation and normalization over the dataset -----------
    ftval = zeros(nt,1);
    for i = 1:nt
        ftval(i) = cases(c).fTraw(Tvec(i));
    end
    ftmax = max(ftval);
    a_term = 1/max(ftmax, eps);      % avoid division by zero
    ftval_norm = a_term .* ftval;   % in [0,1]

    % ---- % Simulate the ODE system ----------------------------------------------
    [E, N, A] = simulate_JL( ...
        Tvec, Hvec, Rvec, Wvec, dates, ht, ...
        ftval_norm, fH, fR, fW, gW, ...
        beta_max, muE_min, alphaT_E, alphaH_E, gammaE_max, ...
        muN_min, kW_mort, gammaN_max, alphaT_N, alphaH_N, ...
        muA_min, betaT_A, betaH_A, ...
        Topt_case, Hopt, Tmin_spost);

    % ----% Save results ---------------------------------------------------
    Results.(thermoName).E = E(:);
    Results.(thermoName).N = N(:);
    Results.(thermoName).A = A(:);
    Results.(thermoName).ftval_norm = ftval_norm(:);
    Results.(thermoName).Topt = Topt_case;

    % ---- Global plot (optional) -------------------------------
    figure
    plot(dates, E, 'r', 'LineWidth', 1.5); hold on
    plot(dates, N, 'g', 'LineWidth', 1.5);
    plot(dates, A, 'b', 'LineWidth', 1.5); hold off
    grid on
    xlabel('Date')
    ylabel('Populations per leaf')
    title(sprintf('J. lybica ODE (%s)', thermoName));
end

%% ===================== OBSERVED DATA 2002 ======================
anno_cmp = 2002;

% Adulti (trap)
obsDatesA = datetime(2002, [5 6 6 7 7 7 8 8 9 9 10 10], ...
                           [22 5 19 3 17 31 14 28 11 25 9 23]);
obsAdultsTrap = [0 0 0 10 35 80 150 140 200 265 230 190];

% Nymphs (leaf, whole community)
obsDatesN = datetime(2002, [5 5 6 6 6 7 7 8 8 9 9 10 10], ...
                           [5 19 2 16 30 14 28 11 25 8 22 6 20]);
obsNymphsAllLeaf = [0.00 0.00 0.08 0.12 0.20 0.25 0.27 0.28 0.61 2.00 4.10 2.00 0.43];

obsDatesA = obsDatesA(:); obsAdultsTrap = obsAdultsTrap(:);
obsDatesN = obsDatesN(:); obsNymphsAllLeaf = obsNymphsAllLeaf(:);

%% ===================== YEAR MASK 2002 (meteo) ===============
maskY = year(dates) == anno_cmp;
datesY = dates(maskY);
TY = Tvec(maskY); HY = Hvec(maskY); WY = Wvec(maskY); RY = Rvec(maskY);

%% ===================== FONT FIGURE ==============================
fsAxes=14; fsTitle=16; fsLabel=15; fsLegend=13;

%% ===================== 2002 COMPARISON FOR BOTH CASES ========
for c = 1:numel(cases)
    thermoName = cases(c).name;
 

    E_full = Results.(thermoName).E;
    N_full = Results.(thermoName).N;
    A_full = Results.(thermoName).A;
    ftval_norm = Results.(thermoName).ftval_norm;

    EY = E_full(maskY); NY = N_full(maskY); AY = A_full(maskY);
    EY = EY(:); NY = NY(:); AY = AY(:); datesY = datesY(:);

    % ---- Interpolate the model onto the observation dates -----------------
    A_mod_obs = interp1(datenum(datesY), AY, datenum(obsDatesA), 'linear', 'extrap');
    N_mod_obs = interp1(datenum(datesY), NY, datenum(obsDatesN), 'linear', 'extrap');
    A_mod_obs = A_mod_obs(:);
    N_mod_obs = N_mod_obs(:);

    % ================== kA (picco vs picco) ======================
    [obsPeakTrap, iObsPeak] = max(obsAdultsTrap);
    tObsPeak = obsDatesA(iObsPeak);

    [modPeakA, iModPeak] = max(AY);
    tModPeak = datesY(iModPeak);

    if isfinite(obsPeakTrap) && obsPeakTrap > 0 && isfinite(modPeakA)
        kA = modPeakA / obsPeakTrap;
    else
        kA = NaN;
        warning('[%s] kA non stimabile.', thermoName);
    end

    obsAdultsTrap_use = obsAdultsTrap;

    if ~isfinite(kA) || kA <= 0
        AY_scaled = nan(size(AY));
        A_mod_obs_scaled = nan(size(A_mod_obs));
    else
        AY_scaled = AY / kA;
        A_mod_obs_scaled = A_mod_obs / kA;
    end

    % ================== kN (picco vs picco) ======================
    [obsPeakN, iObsPeakN] = max(obsNymphsAllLeaf);
    tObsPeakN = obsDatesN(iObsPeakN);

    [modPeakN, iModPeakN] = max(NY);
    tModPeakN = datesY(iModPeakN);

    if isfinite(obsPeakN) && obsPeakN > 0 && isfinite(modPeakN)
        kN = modPeakN / obsPeakN;
    else
        kN = NaN;
        warning('[%s] kN non stimabile.', thermoName);
    end

    obsNymphsLeaf = obsNymphsAllLeaf;

    if ~isfinite(kN) || kN <= 0
        NY_scaled = nan(size(NY));
        N_mod_obs_scaled = nan(size(N_mod_obs));
    else
        NY_scaled = NY / kN;
        N_mod_obs_scaled = N_mod_obs / kN;
    end

    % ================== PEAKS & METRICS ========================
    [NYmax_scaled, iNmax_s] = max(NY_scaled);
    tNmax_scaled = datesY(iNmax_s);

    [AYmax_scaled, iAmax_s] = max(AY_scaled);
    tAmax_scaled = datesY(iAmax_s);

    [obsNmax, iObsNmax] = max(obsNymphsLeaf);
    tObsNmax = obsDatesN(iObsNmax);

    [obsAmaxTrap, iObsAmax] = max(obsAdultsTrap_use);
    tObsAmax = obsDatesA(iObsAmax);

    rmseA = sqrt(mean((A_mod_obs_scaled - obsAdultsTrap_use).^2, 'omitnan'));
    maeA  = mean(abs(A_mod_obs_scaled - obsAdultsTrap_use), 'omitnan');

    rmseN = sqrt(mean((N_mod_obs_scaled - obsNymphsLeaf).^2, 'omitnan'));
    maeN  = mean(abs(N_mod_obs_scaled - obsNymphsLeaf), 'omitnan');

  % ================== TEXT OUTPUT ==========================
fprintf('\n====================================================\n');
fprintf('  Feudo Arancio %d Comparison  |  %s thermal function\n', anno_cmp, thermoName);
fprintf('====================================================\n');

fprintf('\n[%s] --- Estimate kA from peak (amplitude) ---\n', thermoName);
fprintf('Observed peak (trap): %.2f on %s\n', obsPeakTrap, datestr(tObsPeak,'dd-mmm-yyyy'));
fprintf('Model peak A(t) (leaf): %.3f on %s\n', modPeakA, datestr(tModPeak,'dd-mmm-yyyy'));
fprintf('kA = modelPeak/observedPeak = %.4f\n', kA);

fprintf('\n[%s] --- Estimate kN from peak (amplitude) ---\n', thermoName);
fprintf('Observed nymph peak (all leafhoppers, leaf): %.3f on %s\n', obsPeakN, datestr(tObsPeakN,'dd-mmm-yyyy'));
fprintf('Model peak N(t) J. lybica (leaf): %.3f on %s\n', modPeakN, datestr(tModPeakN,'dd-mmm-yyyy'));
fprintf('kN = modelPeak/observedPeak = %.4f\n', kN);

fprintf('\n[%s] J. lybica ADULTS (model scaled -> trap):\n', thermoName);
fprintf('  kA = %.4f | RMSE = %.3f | MAE = %.3f\n', kA, rmseA, maeA);
fprintf('  Model peak (trap, scaled): %.2f on %s\n', AYmax_scaled, datestr(tAmax_scaled,'dd-mmm-yyyy'));
fprintf('  Observed peak (trap):      %.2f on %s\n', obsAmaxTrap, datestr(tObsAmax,'dd-mmm-yyyy'));

fprintf('\n[%s] NYMPHS (observed = community):\n', thermoName);
fprintf('  kN = %.4f | RMSE = %.3f | MAE = %.3f\n', kN, rmseN, maeN);
fprintf('  Model peak (leaf, scaled): %.2f on %s\n', NYmax_scaled, datestr(tNmax_scaled,'dd-mmm-yyyy'));
fprintf('  Observed total peak:       %.2f on %s\n', obsNmax, datestr(tObsNmax,'dd-mmm-yyyy'));

% ================== NYMPHS/ADULTS FIGURES =======================

    figure
    plot(datesY, NY_scaled, 'g', 'LineWidth', 2); hold on
    scatter(obsDatesN, obsNymphsLeaf, 40, 'k', 'filled');
    hold off; grid on
    xlabel('Date','FontSize',fsLabel)
    ylabel('Nymphs per leaf','FontSize',fsLabel)
    title(sprintf('Feudo Arancio 2002 (%s): nymphs', thermoName), 'FontSize', fsTitle)
    legend({'Model (scaled)','Observed'},'Location','best','FontSize',fsLegend)
    set(gca,'FontSize',fsAxes)

    figure
    plot(datesY, AY_scaled, 'b', 'LineWidth', 2); hold on
    scatter(obsDatesA, obsAdultsTrap_use, 60, 'k', 'filled');
    hold off; grid on
    xlabel('Date','FontSize',fsLabel)
    ylabel('Adults per trap','FontSize',fsLabel)
    title(sprintf('Feudo Arancio 2002 (%s): adults', thermoName), 'FontSize', fsTitle)
    legend({'Model (scaled)','Observed (trap)'},'Location','best','FontSize',fsLegend)
    set(gca,'FontSize',fsAxes)

    % ================== EARLY WARNING  ==============
    fprintf('\n=== Early warning | %s | anno %d ===\n', thermoName, anno_cmp);

    early_warning_indices(anno_cmp, ...
        dates, E_full(:), N_full(:), A_full(:), ...
        Tvec, Hvec, Rvec, ...
        ftval_norm(:), fH, fR, gammaE_max, ...
        ht,thermoName);
    
end

%% ===================== ENVIRONMENTAL DRIVERS 2002 ===================
figure
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile
plot(datesY, TY, 'b', 'LineWidth', 1.5); grid on
title('Temperature','FontSize',fsTitle)
ylabel('T (\circC)','FontSize',fsLabel)
set(gca,'FontSize',fsAxes)

nexttile
plot(datesY, HY, 'r', 'LineWidth', 1.5); grid on
title('Relative Humidity','FontSize',fsTitle)
ylabel('H (%)','FontSize',fsLabel)
set(gca,'FontSize',fsAxes)

nexttile
plot(datesY, WY, 'g', 'LineWidth', 1.5); grid on
title('Wind Speed','FontSize',fsTitle)
ylabel('W (m s^{-1})','FontSize',fsLabel)
set(gca,'FontSize',fsAxes)

nexttile
plot(datesY, RY, 'm', 'LineWidth', 1.5); grid on
title('Solar Radiation','FontSize',fsTitle)
ylabel('R (MJ m^{-2} d^{-1})','FontSize',fsLabel)
set(gca,'FontSize',fsAxes)

xlabel(tl,'Date','FontSize',fsLabel)
axEnv = findall(gcf,'type','axes');
linkaxes(axEnv,'x');

%% ==============================================================
% ===================== LOCAL FUNCTIONS ==========================
% ==============================================================

function [E,N,A] = simulate_JL( ...
    Tvec, Hvec, Rvec, Wvec, dates, ht, ...
    ftval_norm, fH, fR, fW, gW, ...
    beta_max, muE_min, alphaT_E, alphaH_E, gammaE_max, ...
    muN_min, kW_mort, gammaN_max, alphaT_N, alphaH_N, ...
    muA_min, betaT_A, betaH_A, ...
    Topt, Hopt, Tmin_spost)

    nt = length(Tvec);

    E = zeros(nt,1);
    N = zeros(nt,1);
    A = zeros(nt,1);

    A_intro = 2/100;
    E_intro = 0;

    adults_introduced_this_year = false;

    for j = 1:nt-1
        T = Tvec(j);
        H = Hvec(j);
        R = Rvec(j);
        W = Wvec(j);

        doy = day(dates(j),'dayofyear');

        % new year: reset
        if doy == 1
            E(j) = 0;
            N(j) = 0;
            A(j) = 0;
            adults_introduced_this_year = false;
        end

        % adult introduction
        if (T > Tmin_spost) && ~adults_introduced_this_year
            A(j) = A(j) + A_intro;
            E(j) = E(j) + E_intro;
            adults_introduced_this_year = true;
        end

        % rates
        beta   = beta_max   * ftval_norm(j) * fH(H) * fR(R);
        gammaE = gammaE_max * ftval_norm(j) * fH(H) * fR(R);

        muE = muE_min + alphaT_E*abs(T - Topt) + alphaH_E*abs(H - Hopt);

        gammaN = gammaN_max * ftval_norm(j) * gW(W);

        muN = muN_min + alphaT_N*abs(T - Topt) + alphaH_N*abs(H - Hopt);
        muN = muN + kW_mort * fW(W);

        muA = muA_min + betaT_A*abs(T - Topt) + betaH_A*abs(H - Hopt);

        % Euler step
        E(j+1) = E(j) + ht * ( beta*A(j)   - (muE + gammaE)*E(j) );
        N(j+1) = N(j) + ht * ( gammaE*E(j) - (muN + gammaN)*N(j) );
        A(j+1) = A(j) + ht * ( gammaN*N(j) - muA*A(j) );
    end
end

% -------- Analytis raw: a*(T-Tmin)*sqrt(Tmax-T) without the T factor
function f = fT_Analytis_raw(T, Tmin, Tmax)
    if T <= Tmin || T >= Tmax
        f = 0.0;
    else
        f = (T - Tmin) * sqrt(max(Tmax - T, 0));
    end
end

% -------- Briere raw: a*T*(T-Tmin)*sqrt(Tmax-T)
function f = fT_Briere_raw(T, Tmin, Tmax)
    if T <= Tmin || T >= Tmax
        f = 0.0;
    else
        f = T * (T - Tmin) * sqrt(max(Tmax - T, 0));
    end
end

% -------- Topt per Briere classica f(T)=T(T-Tmin)sqrt(Tmax-T)
function Topt = briere_Topt_classic(Tmin, Tmax)
    Topt = (4*Tmax + 3*Tmin + sqrt(16*Tmax^2 - 16*Tmax*Tmin + 9*Tmin^2)) / 10;
end



function early_warning_indices(anno_indici, ...
                               dates, E, N, A, ...
                               Tvec, Hvec, Rvec, ...
                               ftval_norm, fH, fR, gammaE_max, ...
                               ht, varargin)

% ---- optional thermal case name (backward compatible)
thermoName = "";
if ~isempty(varargin)
    thermoName = string(varargin{1});
end

if strlength(thermoName) > 0
    caseTag = " (" + thermoName + ")";
else
    caseTag = "";
end

% Computes ENEI, ARI, and CRI for a given year and defines risk based on N(t):
%   - LOW risk   : from t_BASSO = t_EW (first ENEI > ENEI_ew_thr) while N < N_thr_med
%   - MEDIUM risk: N_thr_med <= N < N_thr_alt
%   - HIGH risk  : N >= N_thr_alt
% Risk assessment starts at t_BASSO; the timeline is shown
% (typically) for May--August.

 
%% ---- Year selection -----------------------------------------------
data_inizio = datetime(anno_indici,1,1);
data_fine   = datetime(anno_indici,12,31);

maskYear = (dates >= data_inizio) & (dates <= data_fine);
if ~any(maskYear)
    warning('Year %d is not available in the dataset.', anno_indici);
    return;
end

    idxYear = find(maskYear);
    datesY  = dates(maskYear);
    EY      = E(maskYear);
    NY      = N(maskYear);
    AY      = A(maskYear);

    datesY = datesY(:);
    EY     = EY(:);
    NY     = NY(:);
    AY     = AY(:);

    %% ---- ENEI (eggs -> nymphs), NOT normalized -----------------------
    TY  = Tvec(idxYear);
    HY  = Hvec(idxYear);
    RY  = Rvec(idxYear);
    TY  = TY(:); HY = HY(:); RY = RY(:);

    ftY = ftval_norm(idxYear);
    ftY = ftY(:);

    gammaE_Y = gammaE_max .* ftY .* fH(HY) .* fR(RY);

    Escale = max(EY);
    if Escale <= 0
        Escale = 1;   % numerical safeguard
    end

    ENEI = cumsum( gammaE_Y .* (EY / Escale) ) * ht;   % cumulativo assoluto

    %% ---- ENEI early warning: start of LOW risk ---------------------
    ENEI_ew_thr = 0;                               % ENEI threshold to indicate "activity is present"
    i_EW = find(ENEI > ENEI_ew_thr, 1, 'first');   % first day with ENEI > threshold

    if ~isempty(i_EW)
        t_BASSO = datesY(i_EW);    % start of LOW risk
    else
        t_BASSO = NaT;
    end

    %% ---- Risk thresholds on N(t) ------------------------------------
    N_thr_med = 0.5;   % nymphs/leaf for MEDIUM risk
    N_thr_alt = 2;     % nymphs/leaf for HIGH risk

    if ~isnat(t_BASSO)
        mask_from_BASSO = (datesY >= t_BASSO);
    else
        mask_from_BASSO = false(size(datesY));
    end

    % first date with N >= N_thr_med AFTER t_BASSO
    idx_med = find(mask_from_BASSO & (NY >= N_thr_med), 1, 'first');
    if ~isempty(idx_med)
        t_MED       = datesY(idx_med);
        ENEI_at_med = ENEI(idx_med);
        lead_B_med  = days(t_MED - t_BASSO);
    else
        t_MED       = NaT;
        ENEI_at_med = NaN;
        lead_B_med  = NaN;
    end

    % first date with N >= N_thr_alt AFTER t_BASSO
    idx_alt = find(mask_from_BASSO & (NY >= N_thr_alt), 1, 'first');
    if ~isempty(idx_alt)
        t_ALTO       = datesY(idx_alt);
        ENEI_at_alt  = ENEI(idx_alt);
        lead_B_alt   = days(t_ALTO - t_BASSO);
    else
        t_ALTO       = NaT;
        ENEI_at_alt  = NaN;
        lead_B_alt   = NaN;
    end

  %% ---- Risk index based ONLY on N(t) ------------------------
RiskN      = zeros(size(NY));
riskClassN = strings(size(NY));

for k = 1:length(NY)
    if ~mask_from_BASSO(k)
        RiskN(k)      = 0.0;
        riskClassN(k) = "No risk";
    else
        if NY(k) < N_thr_med
            RiskN(k)      = 0.25;
            riskClassN(k) = "Low risk";
        elseif NY(k) < N_thr_alt
            RiskN(k)      = 0.60;
            riskClassN(k) = "Medium risk";
        else
            RiskN(k)      = 1.00;
            riskClassN(k) = "High risk";
        end
    end
end

   %% ---- Text output ----------------------------------------------
fprintf('\n=== Early warning%s & nymph-based risk, year %d ===\n', caseTag, anno_indici);

% LOW risk
fprintf('\nLOW risk (ENEI > %.2f):\n', ENEI_ew_thr);
if ~isnat(t_BASSO)
    fprintf('  low-risk start (t_BASSO): %s\n', datestr(t_BASSO,'dd-mmm-yyyy'));
else
    fprintf('  ENEI never exceeds the threshold: risk always absent.\n');
end

% MEDIUM risk
fprintf('\nMEDIUM risk: N >= %.2f nymphs/leaf\n', N_thr_med);
if ~isnat(t_MED)
    fprintf('  first date t_MED: %s\n', datestr(t_MED,'dd-mmm-yyyy'));
    fprintf('  ENEI(t_MED) = %.3f\n', ENEI_at_med);
    if ~isnat(t_BASSO)
        fprintf('  lead time t_BASSO -> t_MED: %.1f days\n', lead_B_med);
    end
else
    fprintf('  MEDIUM threshold not reached after t_BASSO.\n');
end

% HIGH risk
fprintf('\nHIGH risk: N >= %.2f nymphs/leaf\n', N_thr_alt);
if ~isnat(t_ALTO)
    fprintf('  first date t_ALTO: %s\n', datestr(t_ALTO,'dd-mmm-yyyy'));
    fprintf('  ENEI(t_ALTO) = %.3f\n', ENEI_at_alt);
    if ~isnat(t_BASSO)
        fprintf('  lead time t_BASSO -> t_ALTO: %.1f days\n', lead_B_alt);
    end
else
    fprintf('  HIGH threshold not reached after t_BASSO.\n');
end

    %% ---- Plot 1: populations per leaf ----------------------------
    figure
    plot(datesY, EY, 'r', 'LineWidth', 2); hold on
    plot(datesY, NY, 'g', 'LineWidth', 2);
    plot(datesY, AY, 'b', 'LineWidth', 2); hold off
    legend('Eggs','Nymphs','Adults','Location','best')
    xlabel('Date')
    ylabel('Individuals per leaf')
   % title(sprintf('J. lybica - %d (populations per leaf)', anno_indici));
    title(sprintf('J. lybica%s - %d (popolazioni per foglia)', caseTag, anno_indici));

    grid on
    ax = gca;
    ax.XTick = data_inizio:calmonths(1):data_fine;
    if isprop(ax,'XAxis') && isprop(ax.XAxis,'TickLabelFormat')
        ax.XAxis.TickLabelFormat = 'MMM-yyyy';
    end
    ax.XTickLabelRotation = 45;

    %% ---- Plot 2: ENEI + date markers -----------------------------
    figure
    plot(datesY, ENEI, 'm', 'LineWidth', 2); hold on
    yline(ENEI_ew_thr, '--k', 'ENEI_{ew}');
    yMax  = max(ENEI);
    yText = 0.95 * yMax;

    % LOW (green)
    if ~isnat(t_BASSO)
        xline(t_BASSO,'--','Color',[0 0.6 0]);
        text(t_BASSO, yText, datestr(t_BASSO,'dd-mmm'), ...
            'Color',[0 0.6 0], 'Rotation',90, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','top', ...
            'FontWeight','bold');
    end

    % MEDIUM (orange)
    if ~isnat(t_MED)
        xline(t_MED,'--','Color',[1 0.5 0]);
        text(t_MED, yText, datestr(t_MED,'dd-mmm'), ...
            'Color',[1 0.5 0], 'Rotation',90, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','top', ...
            'FontWeight','bold');
    end

    % HIGH (red)
    if ~isnat(t_ALTO)
        xline(t_ALTO,'--','Color',[0.8 0 0]);
        text(t_ALTO, yText, datestr(t_ALTO,'dd-mmm'), ...
            'Color',[0.8 0 0], 'Rotation',90, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','top', ...
            'FontWeight','bold');
    end

    hold off
    xlabel('Date')
    ylabel('ENEI (cumulative, absolute)')
    %title(sprintf('ENEI and risk thresholds - %d', anno_indici));
    title(sprintf('ENEI e soglie di rischio%s - %d', caseTag, anno_indici));

    grid on
    ax = gca;
    ax.XTick = data_inizio:calmonths(1):data_fine;
    if isprop(ax,'XAxis') && isprop(ax.XAxis,'TickLabelFormat')
        ax.XAxis.TickLabelFormat = 'MMM-yyyy';
    end
    ax.XTickLabelRotation = 45;
end
