function y = enei_fT(T, p)
% ENEI_FT  Analytis thermal response
%
%     ftval(i)   = (T(i)-Tmin) * sqrt(max(Tmax - T(i), 0))   for Tmin<T<Tmax
%     ftval_norm = ftval / max(max(ftval), eps)


if ~isfield(p, 'fT_norm') || isempty(p.fT_norm)
    error('enei_fT:normMissing', ...
      ['p.fT_norm is not set. Call p = set_fT_normalization(p, Tvec) after ' ...
       'building nominal parameters and after ANY change to T_min or T_max]);
      
end

if T <= p.T_min || T >= p.T_max
    y = 0; return;
end
raw = (T - p.T_min) * sqrt(max(p.T_max - T, 0));
y = raw / p.fT_norm;     
end
