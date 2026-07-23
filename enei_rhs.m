function [dstate, rates] = enei_rhs(state, T, H, R, W, p)
% ENEI_RHS  Single shared right-hand side of the stage-structured model
% (Eq. 1), used by the main simulation, time-resolved sensitivities, endpoint
% elasticities and finite-range robustness analyses. 
%
% INPUTS
%   state : [E; N; A]  (A = TOTAL adult density. Oviposition uses the
%                       effective per-total-adult coefficient beta_max directly.)
%   T,H,R,W : scalar daily forcing values
%   p       : parameter struct (single source; T_opt already derived,
%             H_m shared between f_H and the humidity-stress reference)
%
% OUTPUTS
%   dstate : [dE; dN; dA]  time derivatives (per day)
%   rates  : struct of the effective rates actually used, exposed so the
%            variational/sensitivity system uses  the same numbers.

E = state(1); N = state(2); A = state(3);

% --- neutralization switches (absent => active); only make_temperature_only
%     sets these to 0. Default 1 leaves the full model exactly unchanged. ---
sw_fH = getfield_def(p,'sw_fH',1);
sw_fR = getfield_def(p,'sw_fR',1);
sw_gW = getfield_def(p,'sw_gW',1);
sw_fW = getfield_def(p,'sw_fW',1);

% --- normalized environmental responses ---
fT = enei_fT(T, p);                       % Analytis, normalized to [0,1]
if sw_fH, fH = 1./(1+exp(-p.c_H*(H - p.H_m))); else, fH = 1; end
if sw_fR, fR = 1./(1+exp(-p.c_R*(R - p.R_m))); else, fR = 1; end
if sw_gW, gW = 1./(1 + p.alpha_W*W);           else, gW = 1; end
if sw_fW, fW = W./(W + p.K_W);                 else, fW = 0; end

% --- effective rates (wind made explicit) ---
beta       = p.beta_max   * fT * fH * fR;
gammaE_eff = p.gamma_E_max* fT * fH * fR;

muE        = p.mu_E_min + p.alpha_T*abs(T-p.T_opt) + p.alpha_H*abs(H-p.H_m);
muN_base   = p.mu_N_min + p.delta_T*abs(T-p.T_opt) + p.delta_H*abs(H-p.H_m);
muA        = p.mu_A_min + p.beta_T *abs(T-p.T_opt) + p.beta_H *abs(H-p.H_m);

wind_mortality  = p.k_W * fW;             % explicit wind mortality
muN_eff         = muN_base + wind_mortality;
gammaN_eff      = p.gamma_N_max * fT * gW;% wind slowdown inside development

% --- oviposition: effective per-total-adult coefficient ---
eggProduction = beta * A;

% --- Eq. (1) ---
dE = eggProduction    - (muE + gammaE_eff) * E;
dN = gammaE_eff * E   - (muN_eff + gammaN_eff) * N;
dA = gammaN_eff * N   - muA * A;

dstate = [dE; dN; dA];

rates = struct('beta',beta,'gammaE',gammaE_eff,'gammaN',gammaN_eff, ...
               'muE',muE,'muN',muN_eff,'muA',muA, ...
               'fW',fW,'gW',gW,'wind_mortality',wind_mortality);
end

function v = getfield_def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
