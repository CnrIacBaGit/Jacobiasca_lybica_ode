function E = run_local_elasticities(T, p0, forcing, Tvec)

% RUN_LOCAL_ELASTICITIES  Normalized local elasticities of continuous
% endpoints for every independently-varied parameter.
%
%   Standard positive parameters: central-difference elasticity
%       E_p^Q = [Q(p(1+eps)) - Q(p(1-eps))] / (2*eps*Q(p)),  eps = 0.01
%


eps = 0.01;
hC  = 0.5;   % Celsius half-step for thermal thresholds

endpoints = {'ENEI_1Aug','Nmax'};
thermalSemi = {'T_min','T_max'};
excluded    = {'T_opt','T_move'};

rows = {};
for i = 1:height(T)
    pid = char(T.parameter_id(i));
    if any(strcmp(pid, excluded)), continue; end
    if lower(string(T.independently_varied(i))) ~= "yes", continue; end

    base = p0.(pid);

    if any(strcmp(pid, thermalSemi))
        pPlus  = apply_parameter_scenario(p0, pid, base + hC, Tvec);
        pMinus = apply_parameter_scenario(p0, pid, base - hC, Tvec);
        mode = "semi_per_1C";
        denomStep = 2*hC;    % gives dQ per 1 C after dividing by step
    else
        pPlus  = apply_parameter_scenario(p0, pid, base*(1+eps), Tvec);
        pMinus = apply_parameter_scenario(p0, pid, base*(1-eps), Tvec);
        mode = "elasticity";
        denomStep = [];      % handled per-endpoint below
    end

    Qp = compute_continuous_endpoints(simulate_enei_model(pPlus,  forcing));
    Qm = compute_continuous_endpoints(simulate_enei_model(pMinus, forcing));
    Q0 = compute_continuous_endpoints(simulate_enei_model(p0,     forcing));

    for k = 1:numel(endpoints)
        ep = endpoints{k};
        dQ = Qp.(ep) - Qm.(ep);
        if mode == "elasticity"
            val = dQ / (2*eps * Q0.(ep));                 % dimensionless
        else
            val = dQ / denomStep;                         % dQ per 1 C (absolute)
        end
        rows(end+1,:) = {string(pid), string(ep), mode, val, Q0.(ep)}; %#ok<AGROW>
    end
end

E = cell2table(rows, 'VariableNames', ...
    {'parameter','endpoint','mode','value','Q_nominal'});


E.note = repmat("", height(E), 1);
E.note(E.parameter=="beta_max") = ...
    "effective per-adult rate (beta_eff = per-female fecundity x female fraction)";
end
