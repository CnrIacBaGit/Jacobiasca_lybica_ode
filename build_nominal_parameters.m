function p = build_nominal_parameters(T)
% BUILD_NOMINAL_PARAMETERS  Assemble the nominal parameter struct from the
% validated interval table. 

p = struct();
for i = 1:height(T)
    pid = char(T.parameter_id(i));
    if strcmp(pid,'T_opt'), continue; end   % derived below
    p.(pid) = todouble(T.nominal(i));
end

p.T_opt = (2*p.T_max + p.T_min)/3;

% nymph scaling factor kN = modelled raw peak / observed peak = 658.6032/4.10,
% used ONLY to place model-derived nymph-risk crossings on the observed scale
% (0.5/1.0/2.0 nymphs/leaf) during calibration. It does NOT affect the E,N,A trajectories.
p.kN_scale = 658.6032 / 4.10;
end

function x = todouble(v)
if isnumeric(v), x = double(v); return; end
x = str2double(strtrim(string(v)));
end
