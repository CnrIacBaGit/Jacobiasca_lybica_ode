function p = apply_parameter_scenario(p0, paramName, value, Tvec)
% APPLY_PARAMETER_SCENARIO  Return a copy of p0 with paramName set to value,
% applying all coupling rules encoded in the interval table.
%
% Coupling rules:
%   - T_min or T_max changed -> recompute T_opt = (2*T_max + T_min)/3 AND
%     recompute the observed-max Analytis normalization p.fT_norm over the
%     full record Tvec (the normalization depends on the thresholds and the
%     data, exactly as in the original solver). Tvec MUST be supplied when
%     perturbing a thermal threshold.
%   - H_m changed -> same field used in f_H and the |H-H_m| stress terms, so
%     a single assignment propagates correctly.
%   - T_move changed -> biofix recomputed downstream in the simulator.
%   - T_opt -> never set independently.

if strcmp(paramName,'T_opt')
    error("apply_parameter_scenario:toptIndependent", ...
          "T_opt is derived and must not be set independently.");
end

p = p0;
if ~isfield(p, paramName)
    error("apply_parameter_scenario:unknownParam","Unknown parameter '%s'.", paramName);
end
p.(paramName) = value;

if any(strcmp(paramName, {'T_min','T_max'}))
    p.T_opt = (2*p.T_max + p.T_min)/3;
    if nargin < 4 || isempty(Tvec)
        error("apply_parameter_scenario:needTvec", ...
          ["Perturbing %s changes the observed-max thermal normalization; " ...
           "pass the full-record Tvec so p.fT_norm can be recomputed."], paramName);
    end
    p = set_fT_normalization(p, Tvec);
end

end
