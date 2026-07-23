function P = JL_load_parameters(csvPath)

if nargin < 1 || isempty(csvPath)
    csvPath = 'ENEI_parameter_intervals_updated.csv';
end
if ~isfile(csvPath)
    error('JL_load_parameters:noCSV','Parameter source not found: %s', csvPath);
end

T   = readtable(csvPath, 'TextType','string', 'VariableNamingRule','preserve');
ids = strtrim(string(T.parameter_id));

MAP = { ...
    'beta_max',{'beta_max'}; 'mu_E_min',{'muE_min'}; 'alpha_T',{'alphaT_E'}; ...
    'alpha_H',{'alphaH_E'};  'gamma_E_max',{'gammaE_max'}; 'mu_N_min',{'muN_min'}; ...
    'k_W',{'kW_mort'};       'gamma_N_max',{'gammaN_max'}; 'delta_T',{'alphaT_N'}; ...
    'delta_H',{'alphaH_N'};  'mu_A_min',{'muA_min'};       'beta_T',{'betaT_A'}; ...
    'beta_H',{'betaH_A'};    'T_min',{'Tmin'};             'T_max',{'Tmax'}; ...
    'T_move',{'Tmin_spost'}; 'H_m',{'H0','Hopt'};          'c_H',{'cH'}; ...
    'c_R',{'cR'};            'R_m',{'R0'};                 'K_W',{'KW_half'}; ...
    'alpha_W',{'alphaW_dev'};'A_intro',{'A_intro'} };


P = struct();
for i = 1:size(MAP,1)
    cid = MAP{i,1}; names = MAP{i,2};
    r = find(ids == string(cid));
    if isempty(r)
        error('JL_load_parameters:missingParam','Parameter "%s" not in %s.', cid, csvPath);
    elseif numel(r) > 1
        error('JL_load_parameters:duplicateParam','Parameter "%s" duplicated.', cid);
    end
    v = T.nominal(r);
    if ~isnumeric(v), v = str2double(strtrim(string(v))); else, v = double(v); end
    if isnan(v)
        error('JL_load_parameters:badValue','Parameter "%s" has no numeric nominal.', cid);
    end
    for k = 1:numel(names), P.(names{k}) = v; end
end

P.Topt = (2*P.Tmax + P.Tmin) / 3;    % derived, never read from file

fprintf('[JL_load_parameters] %s: 24 parameters, Topt = %.4f degC.\n', csvPath, P.Topt);
end
