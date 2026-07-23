function p = set_fT_normalization(p, Tvec)
% SET_FT_NORMALIZATION  Compute and store the Analytis normalizing constant


T = Tvec(:);
raw = zeros(numel(T),1);
in = (T > p.T_min) & (T < p.T_max);
raw(in) = (T(in) - p.T_min) .* sqrt(max(p.T_max - T(in), 0));
p.fT_norm = max(max(raw), eps);
end
