function ind = compute_indicator(indicatorName, p, forcingYear)
% COMPUTE_INDICATOR  Return the cumulative indicator trajectory for one year
% on its own absolute scale. 
%
%   'ENEI'      : full multi-driver ENEI (Eq. 8) via simulate_enei_model.
%   'DD'        : cumulative degree-days above T_min, capped at T_max,
%                 

% OUTPUT ind: struct with .dates and .value (cumulative), plus .biofix.

switch indicatorName
    case 'ENEI'
        s = simulate_enei_model(p, forcingYear);
        ind = struct('dates',s.dates,'value',s.ENEI,'biofix',s.biofix);
    case 'DD'
        ind = degree_days(p, forcingYear);
    otherwise
        error('compute_indicator:unknown','Unknown indicator %s.', indicatorName);
end
end

function ind = degree_days(p, F)
T = F.T(:); dates = F.dates(:);
i0 = find(T >= p.T_move, 1, 'first');
if isempty(i0)
    error('degree_days:noBiofix','No biofix (T>=T_move) in year %d.', F.year);
end
n = numel(T);
DD = zeros(n,1);
for j = i0:n
    contrib = max(0, min(T(j), p.T_max) - p.T_min);
    if j == i0, DD(j) = contrib; else, DD(j) = DD(j-1) + contrib; end
end
ind = struct('dates',dates,'value',DD,'biofix',dates(i0));
end
