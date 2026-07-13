function T = load_parameter_intervals(csvPath)
% LOAD_PARAMETER_INTERVALS  Read the single-source interval table and
% validate its invariants. The CSV is the authoritative source; no numeric
% values are duplicated in the scripts.

if nargin < 1
    csvPath = "ENEI_parameter_intervals_updated.csv";
end

T = readtable(csvPath, "TextType","string");

required = ["parameter_id","nominal","lower","upper","independently_varied"];
missing = setdiff(required, string(T.Properties.VariableNames));
if ~isempty(missing)
    error("load_parameter_intervals:missingCols", ...
          "CSV missing required columns: %s", strjoin(missing, ", "));
end

% expected parameter set (T_opt is derived, still listed)
expected = ["T_min","T_max","T_opt","beta_max","mu_E_min","mu_N_min", ...
    "mu_A_min","gamma_E_max","gamma_N_max","A_intro","T_move", ...
    "alpha_T","alpha_H","delta_T","delta_H","beta_T","beta_H", ...
    "K_W","k_W","alpha_W","c_H","H_m","c_R","R_m"];
have = string(T.parameter_id);
missP = setdiff(expected, have);
if ~isempty(missP)
    error("load_parameter_intervals:missingParams", ...
          "CSV missing parameters: %s", strjoin(missP, ", "));
end

% unique ids
if numel(unique(have)) ~= numel(have)
    error("load_parameter_intervals:dupIds","parameter_id values are not unique.");
end

for i = 1:height(T)
    pid = T.parameter_id(i);
    nom = todouble(T.nominal(i));
    lo  = todouble(T.lower(i));
    hi  = todouble(T.upper(i));

    if pid == "T_opt"
        % derived: must have NO bounds and must not be independently varied
        if ~isnan(lo) || ~isnan(hi)
            error("load_parameter_intervals:toptBounds", ...
                  "T_opt must have empty lower/upper (it is derived).");
        end
        if lower(string(T.independently_varied(i))) ~= "no"
            error("load_parameter_intervals:toptVaried", ...
                  "T_opt.independently_varied must be No.");
        end
        continue;
    end

    if isnan(nom)
        error("load_parameter_intervals:nanNominal","%s has non-numeric nominal.", pid);
    end
    % lower <= nominal <= upper where bounds exist
    if ~isnan(lo) && ~isnan(hi)
        if ~(lo <= nom && nom <= hi)
            error("load_parameter_intervals:orderViol", ...
                  "%s violates lower <= nominal <= upper (%.6g, %.6g, %.6g).", ...
                  pid, lo, nom, hi);
        end
    end
end

fprintf("[load_parameter_intervals] %d parameters loaded and validated.\n", height(T));
end

function x = todouble(v)
if isnumeric(v), x = double(v); return; end
s = strtrim(string(v));
if s == "" || ismember(lower(s), ["--","-","na","nan"]) , x = NaN; return; end
x = str2double(s);
end
