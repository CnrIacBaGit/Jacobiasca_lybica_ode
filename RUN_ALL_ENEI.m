% RUN_ALL_ENEI  Master launcher for the complete revised ENEI repository.
%
% Run this file from MATLAB. It changes the current folder to the repository,
% verifies all required project files and inputs, then regenerates:
%   1) main manuscript figures/tables;
%   2) endpoint elasticities and range robustness outputs;
%   3) the  low/nominal/high trajectory figure;
%   4) the five Appendix-A sensitivity figures.
%
% PRINCIPAL MODEL FILE: script_solve_ode_Arancio.m

codeDir = fileparts(mfilename('fullpath'));
cd(codeDir);
addpath(codeDir);

verify_repository();

fprintf('\n[1/4] Main manuscript simulation and calibration...\n');
script_solve_ode_Arancio;

fprintf('\n[2/4] Elasticities, interactions and robustness analyses...\n');
run_enei_analysis('ENEI_parameter_intervals_updated.csv', ...
    'Arancio_dati_2000_2025.mat', 'enei_out');

fprintf('\n[3/4] Reviewer-2 stage-trajectory figure...\n');
make_figure_R22('ENEI_parameter_intervals_updated.csv', ...
    'Arancio_dati_2000_2025.mat', 'figure_R22_stage_trajectories.png');

fprintf('\n[4/4] Appendix-A time-dependent sensitivity figures...\n');
script_sensitivity_aintro();
exportgraphics(gcf, 'arancio_sens_aintro_rev.jpg', 'Resolution', 300);
close(gcf);

script_sensitivity_temperature();
exportgraphics(gcf, 'arancio_sens_tbar_rev.jpg', 'Resolution', 300);
close(gcf);

script_sensitivity_humidity();
exportgraphics(gcf, 'arancio_sens_hbar_rev.jpg', 'Resolution', 300);
close(gcf);

script_sensitivity_radiation();
exportgraphics(gcf, 'arancio_sens_rbar_rev.jpg', 'Resolution', 300);
close(gcf);

script_sensitivity_wind();
exportgraphics(gcf, 'arancio_sens_wbar_rev.jpg', 'Resolution', 300);
close(gcf);

fprintf('\n[RUN_ALL_ENEI] Complete.\n');
fprintf('Main figures: figures_paper/\n');
fprintf('Main tables:  tables_paper/\n');
fprintf('Analysis outputs: enei_out/\n');
