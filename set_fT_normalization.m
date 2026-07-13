function p = set_fT_normalization(p, Tvec)
% SET_FT_NORMALIZATION  Compute and store the Analytis normalizing constant
% EXACTLY as the original solver does: the maximum raw Analytis value over
% the ENTIRE forcing record Tvec (2000-2025), using the current T_min/T_max.
%
%     raw(i)    = (Tvec(i)-T_min)*sqrt(max(T_max - Tvec(i),0))  in-window
%     p.fT_norm = max(max(raw), eps)
%
% MUST be called:
%   - once after build_nominal_parameters, and
%   - again inside apply_parameter_scenario whenever T_min or T_max changes,
%     because the normalization depends on both the data and the thresholds.
%
% Passing the full-record Tvec (not a single year) reproduces the original
% global normalization; the same constant is then used for every year.

T = Tvec(:);
raw = zeros(numel(T),1);
in = (T > p.T_min) & (T < p.T_max);
raw(in) = (T(in) - p.T_min) .* sqrt(max(p.T_max - T(in), 0));
p.fT_norm = max(max(raw), eps);
end
