function make_figure_R22(csvPath, matPath, outFile)
% MAKE_FIGURE_R22  Generate the stage-specific trajectory figure requested by
% Reviewer 2 (comment R2.2): simulated E(t), N(t), A(t) and the resulting ENEI
% signal under the low / nominal / high settings of the parameters the Reviewer
% singled out -- the stage-specific baseline mortalities (Table 2) and the
% maximum egg-hatching and nymph-development rates (Table 3).
%
% The figure is what makes the answer to R2.3 visible: for the mortalities the
% amplitude of the curves changes while the ENEI warning date barely moves,
% whereas for gamma_E,max and gamma_N,max the warning date itself shifts.
%
% USAGE
%   make_figure_R22                                  % uses defaults
%   make_figure_R22(csv, mat, 'fig_R22.png')
%
% Requires the validated pipeline (load_parameter_intervals, load_forcing,
% build_nominal_parameters, set_fT_normalization, apply_parameter_scenario,
% simulate_enei_model, observed_2002).

if nargin<1 || isempty(csvPath), csvPath='ENEI_parameter_intervals_updated.csv'; end
if nargin<2 || isempty(matPath), matPath='Arancio_dati_2000_2025.mat'; end
if nargin<3 || isempty(outFile), outFile='figure_R22_stage_trajectories.png'; end

% ---- setup (same validated pipeline as run_enei_analysis) ----
Tint = load_parameter_intervals(csvPath);
p0   = build_nominal_parameters(Tint);
F    = load_forcing(matPath);
p0   = set_fT_normalization(p0, F.T);
obs  = observed_2002();

idx2002 = year(F.dates)==2002;
F2002 = struct('T',F.T(idx2002),'H',F.H(idx2002),'R',F.R(idx2002), ...
               'W',F.W(idx2002),'dates',F.dates(idx2002),'year',2002);

% ---- the five parameters raised by Reviewer 2 ----
params = {'mu_E_min','mu_N_min','mu_A_min','gamma_E_max','gamma_N_max'};
labels = {'\mu_{E,min}','\mu_{N,min}','\mu_{A,min}','\gamma_{E,max}','\gamma_{N,max}'};

% nominal fixed ENEI threshold (low-risk, L=21) to mark the warning date;
% this is the calibrated operational threshold reported in the manuscript.
theta_low = 18.60;

nP = numel(params);
fig = figure('Position',[100 100 1250 260*nP],'Color','w');

for i = 1:nP
    pid = params{i};
    r   = find(Tint.parameter_id == string(pid), 1);
    lo  = str2double(string(Tint.lower(r)));
    hi  = str2double(string(Tint.upper(r)));
    nomv = p0.(pid);
    vals = [lo nomv hi];
    styles = {'--','-',':'};
    cols   = [0.20 0.40 0.80; 0.10 0.10 0.10; 0.85 0.30 0.15];
    scnLab = {'low','nominal','high'};

    S = cell(1,3);
    for s = 1:3
        p = apply_parameter_scenario(p0, pid, vals(s), F.T);
        S{s} = simulate_enei_model(p, F2002);
    end

    % --- E, N, A panels + ENEI panel ---
    names = {'E','N','A'};
    ylabs = {'Eggs (per leaf)','Nymphs (per leaf)','Adults (per leaf)'};
    for k = 1:3
        subplot(nP,4,(i-1)*4+k); hold on; box on
        for s = 1:3
            plot(S{s}.dates, S{s}.(names{k}), styles{s}, ...
                 'Color', cols(s,:), 'LineWidth', 1.4);
        end
        xlim([datetime(2002,6,1) datetime(2002,10,15)]);
        ylabel(ylabs{k},'FontSize',8);
        title(sprintf('%s  |  %s', labels{i}, names{k}),'FontSize',9);
        if i==nP, xlabel('2002','FontSize',8); end
        set(gca,'FontSize',8);
    end

    % --- ENEI + warning dates: this is the panel that answers R2.3 ---
    subplot(nP,4,(i-1)*4+4); hold on; box on
    for s = 1:3
        plot(S{s}.dates, S{s}.ENEI, styles{s}, 'Color', cols(s,:), 'LineWidth',1.4);
        kx = find(S{s}.ENEI >= theta_low, 1, 'first');
        if ~isempty(kx)
            xline(S{s}.dates(kx), styles{s}, 'Color', cols(s,:), ...
                  'LineWidth', 1.0, 'Alpha', 0.6);
        end
    end
    yline(theta_low,'k-','LineWidth',1.0,'Alpha',0.5);
    xlim([datetime(2002,6,1) datetime(2002,10,15)]);
    ylim([0 400]);
    ylabel('ENEI','FontSize',8);
    title(sprintf('%s  |  ENEI + warning date', labels{i}),'FontSize',9);
    if i==nP, xlabel('2002','FontSize',8); end
    if i==1
        legend(scnLab,'Location','northwest','FontSize',7,'Box','off');
    end
    set(gca,'FontSize',8);
end

sgtitle({'Stage-specific trajectories under low / nominal / high parameterizations (2002)', ...
         'Vertical lines: ENEI crossing of the calibrated low-risk warning threshold'}, ...
         'FontSize',11,'FontWeight','bold');

exportgraphics(fig, outFile, 'Resolution', 300);
fprintf('[make_figure_R22] written: %s\n', outFile);

% ---- also print the numbers the figure illustrates ----
fprintf('\nWarning-date shift at the low-risk threshold (theta = %.2f):\n', theta_low);
fprintf('%-14s %-12s %-12s %-12s\n','parameter','low','nominal','high');
for i = 1:nP
    pid = params{i};
    r   = find(Tint.parameter_id == string(pid), 1);
    vals = [str2double(string(Tint.lower(r))), p0.(pid), str2double(string(Tint.upper(r)))];
    dts = strings(1,3);
    for s = 1:3
        p = apply_parameter_scenario(p0, pid, vals(s), F.T);
        S1 = simulate_enei_model(p, F2002);
        kx = find(S1.ENEI >= theta_low, 1, 'first');
        if isempty(kx), dts(s) = "not reached";
        else, dts(s) = string(datestr(S1.dates(kx),'dd-mmm')); end
    end
    fprintf('%-14s %-12s %-12s %-12s\n', pid, dts(1), dts(2), dts(3));
end
end
