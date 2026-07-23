function p = apply_parameter_scenario(p0, paramName, value, Tvec)
% APPLY_PARAMETER_SCENARIO  Return a copy of p0 with paramName set to value,
% applying all coupling rules encoded in the interval table.
%


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
