function run_enei_analysis(csvPath, matPath, outDir)
% RUN_ENEI_ANALYSIS  End-to-end driver on the real 2000-2025 forcing.
% Order (agreed):
%   load/clean forcing -> annual simulation -> rolling calibration
%   -> 2002 evaluation -> range robustness. Temperature-only runs in the
%   SAME cycle with its OWN absolute thresholds and identical protocol.
%
if nargin < 3, outDir = "enei_out"; end
if nargin < 2, matPath = "Arancio_dati_2000_2025.mat"; end
if nargin < 1, csvPath = "ENEI_parameter_intervals_updated.csv"; end

Tint = load_parameter_intervals(csvPath);
p0   = build_nominal_parameters(Tint);
F    = load_forcing(matPath);
obs  = observed_2002();

% observed-max Analytis normalization over the FULL record, as in the
% original solver. Recomputed inside apply_parameter_scenario when a thermal
% threshold is perturbed (Tvec passed through there).
p0   = set_fT_normalization(p0, F.T);

% nominal preflight on 2002 slice
F2002 = struct('T',F.T(year(F.dates)==2002),'H',F.H(year(F.dates)==2002), ...
               'R',F.R(year(F.dates)==2002),'W',F.W(year(F.dates)==2002), ...
               'dates',F.dates(year(F.dates)==2002),'year',2002);
preflight_checks(p0, F2002);

% continuous-endpoint elasticities (finite differences via simulate; no Jacobian)
elast = run_local_elasticities(Tint, p0, F2002, F.T);

% interaction checks
interactions = run_interaction_checks(p0, F2002);

% range robustness. PRIMARY = Protocol B (fixed nominal thresholds), which
% isolates warning-timing sensitivity from recalibration noise. Protocol A
% (end-to-end recalibration) is kept as a calibration-stability diagnostic.
indicators = {'ENEI','DD'};
rrB = struct();   % Protocol B, primary
rrA = struct();   % Protocol A, diagnostic
for k = 1:numel(indicators)
    fld = matlab.lang.makeValidName(indicators{k});
    rrB.(fld) = run_range_robustness_B(Tint, p0, F, obs, 21, 2, indicators{k});
    rrA.(fld) = run_range_robustness(Tint, p0, F, obs, 21, 2, indicators{k});
end
rr = struct('protocolB_primary', rrB, 'protocolA_diagnostic', rrA);

export_enei_results(outDir, elast, interactions, rr);
fprintf('[run_enei_analysis] complete. Elasticities, interactions, and\n');
fprintf('  range robustness for ENEI / DD written to %s.\n', outDir);
end
