function Q = compute_continuous_endpoints(sim)
% COMPUTE_CONTINUOUS_ENDPOINTS  Differentiable endpoints used for elasticity.
%   Q.ENEI_1Aug : cumulative ENEI at the FIXED calendar date 1 August 2002.
%                 This date coincides with the nominal low-risk warning phase
%                 but is held unchanged under all parameter perturbations.
%   Q.Nmax      : seasonal peak nymph density.
%
% Timing endpoints (nymph-peak date, warning dates, lead times) are NOT here;
% they are handled in days by compute_warning_metrics, because they are
% threshold-crossing / argmax quantities and not smoothly differentiable.

yr = year(sim.dates(1));
tstar = datetime(yr, 8, 1);
i = find(sim.dates <= tstar, 1, 'last');
if isempty(i), i = 1; end

Q = struct();
Q.ENEI_1Aug = sim.ENEI(i);
Q.Nmax = sim.Nmax;
end
