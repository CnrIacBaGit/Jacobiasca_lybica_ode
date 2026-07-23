function script_sensitivity_aintro(csvPath, matPath)
% SCRIPT_SENSITIVITY_AINTRO  sensitivity to the initial adult pulse (A_{\\mathrm{intro}}).
%


if nargin < 1 || isempty(csvPath), csvPath = 'ENEI_parameter_intervals_updated.csv'; end
if nargin < 2 || isempty(matPath), matPath = 'Arancio_dati_2000_2025.mat'; end

Tint = load_parameter_intervals(csvPath);
p    = build_nominal_parameters(Tint);
F    = load_forcing(matPath);
p    = set_fT_normalization(p, F.T);      % observed-max over the FULL record

y    = year(F.dates) == 2002;
F2002 = struct('T',F.T(y),'H',F.H(y),'R',F.R(y),'W',F.W(y), ...
               'dates',F.dates(y),'year',2002);

S = enei_sensitivity('aintro', p, F2002);

figure
plot(S.dates, S.sE, 'r', 'LineWidth', 2); hold on
plot(S.dates, S.sN, 'g', 'LineWidth', 2)
plot(S.dates, S.sA, 'b', 'LineWidth', 2)
ylabel('Sensitivity')
title('Sensitivity to $A_{\\mathrm{intro}}$', 'Interpreter', 'latex');
grid on
legend('s_{E}','s_{N}','s_{A}', 'Location', 'best')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

fprintf('[aintro] peak sE = %.4g | sN = %.4g | sA = %.4g\n', ...
    S.sE(end), S.sN(end), S.sA(end));
end
