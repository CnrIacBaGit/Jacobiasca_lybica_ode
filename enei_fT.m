function y = enei_fT(T, p)
% ENEI_FT  Analytis thermal response, normalized EXACTLY as in the original
% solver (script_solve_ode_Arancio.m, lines 235-239):
%
%     ftval(i)   = (T(i)-Tmin) * sqrt(max(Tmax - T(i), 0))   for Tmin<T<Tmax
%     ftval_norm = ftval / max(max(ftval), eps)
%
% i.e. the raw Analytis curve is divided by the MAXIMUM RAW VALUE OBSERVED
% over the whole forcing record, NOT by the theoretical peak at T_opt. The
% two differ whenever the observed temperatures do not reach exactly
% T_opt = (2*Tmax+Tmin)/3, and the observed-max normalization can exceed 1.
%
% Because the normalizing constant depends on the data AND on Tmin/Tmax, it
% must be recomputed whenever Tmin or Tmax is perturbed. It is stored in
% p.fT_norm by set_fT_normalization(p, Tvec); this function requires that
% field and errors if it is missing, to prevent silent use of a wrong
% (e.g. theoretical) normalization.

if ~isfield(p, 'fT_norm') || isempty(p.fT_norm)
    error('enei_fT:normMissing', ...
      ['p.fT_norm is not set. Call p = set_fT_normalization(p, Tvec) after ' ...
       'building nominal parameters and after ANY change to T_min or T_max, ' ...
       'so the observed-max normalization matches the original solver.']);
end

if T <= p.T_min || T >= p.T_max
    y = 0; return;
end
raw = (T - p.T_min) * sqrt(max(p.T_max - T, 0));
y = raw / p.fT_norm;     % NOT clipped to 1: original does not clip either
end
