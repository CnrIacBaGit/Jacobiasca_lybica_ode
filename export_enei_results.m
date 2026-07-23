function export_enei_results(outDir, elast, interactions, rangeRob)
% EXPORT_ENEI_RESULTS  Write analysis outputs to disk.


if ~exist(outDir,'dir'), mkdir(outDir); end

if nargin >= 2 && ~isempty(elast)
    writetable(elast, fullfile(outDir,'local_elasticities.csv'));
end

if nargin >= 3 && ~isempty(interactions)
    % effective-flux range is the headline number for the confounding discussion
    fid = fopen(fullfile(outDir,'interaction_summary.txt'),'w');
    fprintf(fid, "beta_eff x A_intro product range: [%.4g, %.4g]\n", ...
        interactions.grid_product_range(1), interactions.grid_product_range(2));
    fclose(fid);
    save(fullfile(outDir,'interactions.mat'), '-struct', 'interactions');
end

if nargin >= 4 && ~isempty(rangeRob)
    save(fullfile(outDir,'range_robustness.mat'), '-struct', 'rangeRob');

  
    write_protocol_csv(outDir, 'range_robustness_protocolB.csv', ...
                       rangeRob, 'protocolB_primary');
    write_protocol_csv(outDir, 'range_robustness_protocolA.csv', ...
                       rangeRob, 'protocolA_diagnostic');
end

fprintf("[export_enei_results] written to %s\n", outDir);
end

% ---------------------------------------------------------------------------
function write_protocol_csv(outDir, fname, rangeRob, protocolField)
% Flatten the per-indicator tables (ENEI, DD) of one protocol into a single
% CSV, tagging each row with its indicator so the two can be told apart.
if ~isfield(rangeRob, protocolField), return; end
P = rangeRob.(protocolField);

out = table();
inds = fieldnames(P);
for k = 1:numel(inds)
    Tk = P.(inds{k});
    if isempty(Tk), continue; end
    Tk.indicator_group = repmat(string(inds{k}), height(Tk), 1);
    if isempty(out)
        out = Tk;
    else
        out = [out; Tk]; %#ok<AGROW>
    end
end
if isempty(out), return; end

writetable(out, fullfile(outDir, fname));
fprintf("[export_enei_results]   %s (%d rows)\n", fname, height(out));
end
