function preflight_checks(p, forcing)
% PREFLIGHT_CHECKS  Hard gates that must pass before any analysis runs.
% Any violation raises an error and stops execution, so the analyses can
% never be built on an inconsistent model implementation.

tol = 1e-9;

% (1) wind mortality term must be present and active in the nymph equation.
%     Probe the simulator with two wind levels; N-dynamics must differ.
f1 = forcing; f2 = forcing;
f2.W = forcing.W * 3 + 2;                 % markedly higher wind
o1 = simulate_enei_model(p, f1);
o2 = simulate_enei_model(p, f2);
if abs(o1.Nmax - o2.Nmax) < 1e-6
    error("preflight:windTermMissing", ...
      ["Nymph peak is insensitive to wind: the k_W*f_W(W) mortality term " ...
       "(and/or g_W development slowdown) appears to be missing from the " ...
       "nymph update."]);
end

% (2) T_opt must equal the derived formula (no stale independent value).
Topt_expected = (2*p.T_max + p.T_min)/3;
if abs(p.T_opt - Topt_expected) > tol
    error("preflight:toptStale", ...
      "T_opt=%.6g but (2*T_max+T_min)/3=%.6g; recompute after perturbation.", ...
      p.T_opt, Topt_expected);
end

% (3) no negative states beyond numerical tolerance in the nominal run.
negtol = -1e-8;
if any(o1.E < negtol) || any(o1.N < negtol) || any(o1.A < negtol)
    error("preflight:negativeStates", ...
      "Negative state values beyond tolerance; check step size / rates.");
end


fprintf("[preflight_checks] all gates passed.\n");
end
