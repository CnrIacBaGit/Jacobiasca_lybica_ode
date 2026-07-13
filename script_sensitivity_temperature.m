function script_sensitivity_temperature(csvPath, matPath)
% SCRIPT_SENSITIVITY_TEMPERATURE  sensitivity to the mean temperature (\\bar{T}).
%
% This script no longer contains the model equations or any parameter value.
% The state is advanced by ENEI_RHS, the shared routine used by the
% sensitivity, elasticity and robustness analyses, and the variational system
% is built on the rates ENEI_RHS returns (see enei_sensitivity.m). The main
% simulation reads from the same parameter table and uses algebraically
% identical rates. The wind contribution is visible immediately before the
% state update (r.muN contains k_W*f_W; r.gammaN contains g_W).

if nargin < 1 || isempty(csvPath), csvPath = 'ENEI_parameter_intervals_updated.csv'; end
if nargin < 2 || isempty(matPath), matPath = 'Arancio_dati_2000_2025.mat'; end

Tint = load_parameter_intervals(csvPath);
p    = build_nominal_parameters(Tint);
F    = load_forcing(matPath);
p    = set_fT_normalization(p, F.T);      % observed-max over the FULL record

y    = year(F.dates) == 2002;
F2002 = struct('T',F.T(y),'H',F.H(y),'R',F.R(y),'W',F.W(y), ...
               'dates',F.dates(y),'year',2002);

S = enei_sensitivity('T', p, F2002);

figure
plot(S.dates, S.sE, 'r', 'LineWidth', 2); hold on
plot(S.dates, S.sN, 'g', 'LineWidth', 2)
plot(S.dates, S.sA, 'b', 'LineWidth', 2)
ylabel('Sensitivity')
title('Sensitivity to $\\bar{T}$', 'Interpreter', 'latex');
grid on
legend('s_{E}','s_{N}','s_{A}', 'Location', 'best')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

fprintf('[temperature] peak sE = %.4g | sN = %.4g | sA = %.4g\n', ...
    S.sE(end), S.sN(end), S.sA(end));
end
