function out = simulate_enei_model(p, forcing, opts)
% SIMULATE_ENEI_MODEL  Integrate the stage-structured model via enei_rhs.
% Within the auxiliary analysis pipeline, the equations live in enei_rhs and
% this routine only marches Euler and accumulates ENEI. Temperature-only runs
% are obtained by passing a struct
% built with make_temperature_only(p) as p -- no branches here, no second
% copy of the model.
%
% INPUTS
%   p       : parameter struct (single source). A = TOTAL adults; oviposition
%             uses the effective per-total-adult coefficient beta_eff directly.
%   forcing : struct .T .H .R .W (daily), .dates (datetime), .year
%   opts    : optional .t0_override (datetime) to force the biofix date
%
% OUTPUT struct: .dates .E .N .A .gammaE .ENEI .biofix .Nmax .tNpeak

if nargin < 3, opts = struct(); end

T = forcing.T(:); H = forcing.H(:); R = forcing.R(:); W = forcing.W(:);
dates = forcing.dates(:);
nDays = numel(T);

% --- biofix: first day T >= T_move, or override ---
if isfield(opts,'t0_override') && ~isempty(opts.t0_override)
    i0 = find(dates >= opts.t0_override, 1, 'first');
else
    i0 = find(T >= p.T_move, 1, 'first');
end
if isempty(i0)
    error('simulate_enei_model:noBiofix', ...
          'No day satisfies T >= T_move = %.3g for year %d.', p.T_move, forcing.year);
end
biofix = dates(i0);

E = zeros(nDays,1); N = zeros(nDays,1); A = zeros(nDays,1);
gammaE_series = zeros(nDays,1);
ENEI = zeros(nDays,1);

A(i0) = p.A_intro;                 % immigration pulse (total adults)

for j = i0:nDays-1
    [dstate, rates] = enei_rhs([E(j); N(j); A(j)], T(j), H(j), R(j), W(j), p);
    E(j+1) = E(j) + dstate(1);     % explicit Euler, unit daily step
    N(j+1) = N(j) + dstate(2);
    A(j+1) = A(j) + dstate(3);

    gammaE_series(j) = rates.gammaE;
end

% ENEI(t) = cumulative egg-to-nymph flux (Eq. 8), computed EXACTLY as the
% original solver: ENEI = cumsum(gammaE .* E) * ht, with ht = 1 (daily step).
% This accumulates the day-j flux AT day j (no one-day index offset), matching
% script_solve_ode_Arancio.m line 897.
ENEI = cumsum(gammaE_series .* E);

[Nmax, iPk] = max(N);

out = struct('dates',dates,'E',E,'N',N,'A',A, ...
             'gammaE',gammaE_series,'ENEI',ENEI, ...
             'biofix',biofix,'Nmax',Nmax,'tNpeak',dates(iPk));
end
