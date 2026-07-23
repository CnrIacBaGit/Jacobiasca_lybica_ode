function S = enei_sensitivity(mode, p, forcing)
% ENEI_SENSITIVITY  Differential (variational) sensitivity of the stage-
% structured model with respect to the initial adult pulse A_intro or to one
% of the mean environmental forcings (Tbar, Hbar, Rbar, Wbar).
%
%   S = enei_sensitivity(mode, p, forcing)
%   mode : 'aintro' | 'T' | 'H' | 'R' | 'W'
%   p    : parameter struct from build_nominal_parameters + set_fT_normalization
%   forcing : struct .T .H .R .W .dates (one year)
%


mode = validatestring(mode, {'aintro','T','H','R','W'});

T = forcing.T(:); H = forcing.H(:); R = forcing.R(:); W = forcing.W(:);
dates = forcing.dates(:);


i0 = find(T >= p.T_move, 1, 'first');
if isempty(i0)
    error('enei_sensitivity:noBiofix','No day with T >= T_move = %.3g.', p.T_move);
end
iEnd = find(T >= p.T_move, 1, 'last');
idx  = (i0:iEnd).';
n    = numel(idx);

% mean forcings over the active window (the "parameter" being perturbed)
Tbar = mean(T(idx)); Hbar = mean(H(idx));
Rbar = mean(R(idx)); Wbar = mean(W(idx));

E = zeros(n,1); N = zeros(n,1); A = zeros(n,1);
sE = zeros(n,1); sN = zeros(n,1); sA = zeros(n,1);

A(1) = p.A_intro;                       % immigration pulse at the biofix
if strcmp(mode,'aintro')
    sA(1) = 1;                          % dA(t0)/dA_intro = 1
end

for j = 1:n-1
    k = idx(j);
    % the perturbed driver is held at its mean; the others keep their daily values
    Tj = T(k); Hj = H(k); Rj = R(k); Wj = W(k);
    switch mode
        case 'T', Tj = Tbar;
        case 'H', Hj = Hbar;
        case 'R', Rj = Rbar;
        case 'W', Wj = Wbar;
    end

    % ---- STATE: advanced by the shared model routine ----
    [dstate, r] = enei_rhs([E(j); N(j); A(j)], Tj, Hj, Rj, Wj, p);
    E(j+1) = E(j) + dstate(1);
    N(j+1) = N(j) + dstate(2);
    A(j+1) = A(j) + dstate(3);

    % ---- parameter derivatives of the SAME rates ----
    d = rate_derivatives(mode, Tj, Hj, Rj, Wj, p, r);

    % ---- variational system, Eq. (A.2)-(A.5) ----
    % note: r.gammaN already contains g_W (Eq. 6), so it enters LINEARLY here.
    sE(j+1) = sE(j) + ( d.beta*A(j) - (d.muE + d.gammaE)*E(j) ...
                        + r.beta*sA(j) - (r.muE + r.gammaE)*sE(j) );
    sN(j+1) = sN(j) + ( d.gammaE*E(j) - (d.muN + d.gammaN)*N(j) ...
                        + r.gammaE*sE(j) - (r.muN + r.gammaN)*sN(j) );
    sA(j+1) = sA(j) + ( d.gammaN*N(j) - d.muA*A(j) ...
                        + r.gammaN*sN(j) - r.muA*sA(j) );
end

S = struct('mode',mode,'dates',dates(idx),'E',E,'N',N,'A',A, ...
           'sE',sE,'sN',sN,'sA',sA, ...
           'Tbar',Tbar,'Hbar',Hbar,'Rbar',Rbar,'Wbar',Wbar, ...
           'biofix',dates(i0));
end

% ---------------------------------------------------------------------------
function d = rate_derivatives(mode, T, H, R, W, p, r)
% Analytical d(rate)/dc for c = A_intro, Tbar, Hbar, Rbar or Wbar.
% Built on the response functions of enei_rhs so the two cannot drift apart.

d = struct('beta',0,'gammaE',0,'gammaN',0,'muE',0,'muN',0,'muA',0);
if strcmp(mode,'aintro')
    return                              % A_intro enters only the initial condition
end

fT = enei_fT(T, p);
fH = 1./(1+exp(-p.c_H*(H - p.H_m)));
fR = 1./(1+exp(-p.c_R*(R - p.R_m)));
gW = 1./(1 + p.alpha_W*W);

switch mode
    case 'T'
        % d/dT of the OBSERVED-MAX normalized Analytis response
        if T > p.T_min && T < p.T_max
            s   = sqrt(max(p.T_max - T, 0));
            dfT = ( s - (T - p.T_min)/(2*s) ) / p.fT_norm;
        else
            dfT = 0;
        end
        d.beta   = p.beta_max    * dfT * fH * fR;
        d.gammaE = p.gamma_E_max * dfT * fH * fR;
        d.gammaN = p.gamma_N_max * dfT * gW;
        sg = sign(T - p.T_opt);
        d.muE = p.alpha_T*sg;  d.muN = p.delta_T*sg;  d.muA = p.beta_T*sg;

    case 'H'
        dfH = p.c_H * exp(-p.c_H*(H - p.H_m)) * fH^2;
        d.beta   = p.beta_max    * fT * dfH * fR;
        d.gammaE = p.gamma_E_max * fT * dfH * fR;
        sg = sign(H - p.H_m);
        d.muE = p.alpha_H*sg;  d.muN = p.delta_H*sg;  d.muA = p.beta_H*sg;

    case 'R'
        dfR = p.c_R * exp(-p.c_R*(R - p.R_m)) * fR^2;
        d.beta   = p.beta_max    * fT * fH * dfR;
        d.gammaE = p.gamma_E_max * fT * fH * dfR;

    case 'W'
        dgW = -p.alpha_W / (1 + p.alpha_W*W)^2;
        dfW =  p.K_W / (W + p.K_W)^2;
        d.gammaN = p.gamma_N_max * fT * dgW;   % wind slows development
        d.muN    = p.k_W * dfW;                % wind adds mortality
end
end
