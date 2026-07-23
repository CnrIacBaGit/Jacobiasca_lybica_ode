function Q = compute_continuous_endpoints(sim)
% COMPUTE_CONTINUOUS_ENDPOINTS  Differentiable endpoints used for elasticity.
%   Q.ENEI_1Aug : cumulative ENEI at the FIXED calendar date 1 August 2002.
%
%   Q.Nmax      : seasonal peak nymph density.
%


yr = year(sim.dates(1));
tstar = datetime(yr, 8, 1);
i = find(sim.dates <= tstar, 1, 'last');
if isempty(i), i = 1; end

Q = struct();
Q.ENEI_1Aug = sim.ENEI(i);
Q.Nmax = sim.Nmax;
end
