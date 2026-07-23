function IC = run_interaction_checks(p0, forcing)
% RUN_INTERACTION_CHECKS 


IC = struct();

% ---- (1) beta_eff x A_intro grid ----
be_grid  = linspace(0.6, 1.4, 5) * p0.beta_max;   % +/-40% around nominal beta_eff
ain_grid = linspace(0.01, 0.04, 5);
[BE, AIN] = meshgrid(be_grid, ain_grid);
Nmax = zeros(size(BE)); ENEI1 = zeros(size(BE)); prod_ = BE.*AIN;
for a = 1:numel(BE)
    p = p0; p.beta_max = BE(a); p.A_intro = AIN(a);
    s = simulate_enei_model(p, forcing);
    Q = compute_continuous_endpoints(s);
    Nmax(a) = Q.Nmax; ENEI1(a) = Q.ENEI_1Aug;
end
IC.grid_betaEff = BE; IC.grid_Aintro = AIN; IC.grid_product = prod_;
IC.grid_Nmax = Nmax; IC.grid_ENEI1Aug = ENEI1;
IC.grid_product_range = [min(prod_(:)) max(prod_(:))];

% ---- (2) T_move x A_intro grid ----
tmove_grid = [20 22 24];
ain3 = [0.01 0.02 0.04];
[TM, A3] = meshgrid(tmove_grid, ain3);
biofix = NaT(size(TM)); gN = zeros(size(TM));
for a = 1:numel(TM)
    p = p0; p.T_move = TM(a); p.A_intro = A3(a);
    s = simulate_enei_model(p, forcing);
    biofix(a) = s.biofix; gN(a) = s.Nmax;
end
IC.tmove_grid = TM; IC.tmove_Aintro = A3;
IC.tmove_biofix = biofix; IC.tmove_Nmax = gN;
end
