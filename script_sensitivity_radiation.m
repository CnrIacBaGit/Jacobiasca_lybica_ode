% Jacobiasca lybica Stage-Structured Model
% Sensitivity Analysis to solar radiation Rbar

clc
clear all
close all

%% Data temperature, relative humidity, solar radiation, wind
% Tvec, Hvec, Rvec, Wvec
load('Arancio_dati_2000_2025.mat')

%% Time interval
year = 2002; % Reference year
startDate = datetime(year,1,1); % First day in the datase
endDate = datetime(year,12,31);
dates     = [startDate:endDate];  
ndate = length(dates);
% Time step
nday = 1;
ht = 1/nday;  

%% Biological parameters
% Reproduction and eggs
beta_max   = 5.0;      % [eggs d^-1 adult^-1]
muE_min    = 0.05;     % [d^-1]
alphaT_E   = 0.02;     % [d^-1 °C^-1]
alphaH_E   = 0.003;    % [d^-1 %^-1]
gammaE_max = 0.22;     % [d^-1]
% Ninfe
muN_min    = 0.04;     % [d^-1]
kW_mort    = 0.10;     % [d^-1]
gammaN_max = 0.10;     % [d^-1]
alphaT_N   = 0.006;    % [d^-1 °C^-1]
alphaH_N   = 0.0012;   % [d^-1 %^-1]
% Adults
muA_min  = 0.03;       % [d^-1]
betaT_A  = 0.015;      % [d^-1 °C^-1]
betaH_A  = 0.002;      % [d^-1 %^-1] 

%% Environmental factors
% Temperature
Tmin       = 11.0;     % [°C] 
Tmax       = 35.0;     % [°C] 
Tmin_spost = 22;       % [°C] Tmove
Topt = (4*Tmax + 3*Tmin + sqrt(16*Tmax^2 - 16*Tmax*Tmin + 9*Tmin^2))/10; % Briere 
% Humidity
cH   = 0.08;           % [adim per %] 
H0   = 60.0;           % [%] 
Hopt = 60.0;           % [%] 
% Solar radiatione
cR = 0.25;             % [adim per MJ m^-2 d^-1] 
R0 = 18.0;             % [MJ m^-2 d^-1] 
% Wind
KW_half    = 3.0;      % [m s^-1]
alphaW_dev = 0.10;     % [adim]

%% Response function
% f_T(T): Briere normalized with respect to the maximum value in the reference year
baseDate  = datetime(2000,1,1);    
offsetDays = days(startDate - baseDate);
idxStart = offsetDays + 1;
idxEnd   = idxStart + ndate - 1;
% Extract values in the reference year
T_year = Tvec(idxStart:idxEnd);
H_year = Hvec(idxStart:idxEnd);
R_year = Rvec(idxStart:idxEnd);
W_year = Wvec(idxStart:idxEnd);
trova = find(T_year > Tmin_spost);
idx1 = trova(1);
idx2 = trova(end);
trova = [idx1:idx2];
ngiorni = length(trova);
nt = ngiorni/ht;

ftval = zeros(ngiorni,1);
for i = 1:ngiorni
    val = trova(i);
    T = T_year(val);
    ftval(i) = fT(T,Tmin,Tmax);
end
ftmax      = max(ftval);
a_briere   = 1/ftmax;                  
ftval_norm = a_briere .* ftval;        % normalization

% f_H(H)
fH = @(H) 1.0 ./ (1.0 + exp(-cH*(H - H0)));

% f_R(R)
fR = @(R) 1.0 ./ (1.0 + exp(-cR*(R - R0)));
dfR = @(R) cR.*exp(-cR*(R - R0)).*fR(R).^2;

% f_W(W) 
fW = @(W) W ./ (W + KW_half);

% g_W(W) 
gW = @(W) 1 ./ (1 + alphaW_dev*W);

%% Inizialization
E = zeros(nt,1);   % eggs
N = zeros(nt,1);   % nymphs 
A = zeros(nt,1);   % adults
sER = zeros(nt,1); % sensitivity eggs
sNR = zeros(nt,1); % sensitivity nymphs
sAR = zeros(nt,1); % sensitivity adults

A_intro       = 2/100;          % adults immigration

adults_introduced_this_year = false;

dates = dates(trova);

% Fixed parameter Rbar
R = mean(R_year(trova));
dfrval = dfR(R);
frval = fR(R);
%% Solve ODE (explicit Euler)
for j = 1:nt-1
    jj = ceil(j*ht);
    doy = day(dates(jj),'dayofyear');
    if doy == 1
        E(j) = 0;
        N(j) = 0;
        A(j) = 0;
        adults_introduced_this_year = false;
    end
    % adults immigration
    if (T >= Tmin_spost) && ~adults_introduced_this_year
        A(j) = A(j) + A_intro;
        adults_introduced_this_year = true;
    end
    valore = trova(j);
    T = T_year(valore);
    H = H_year(valore);
    W = W_year(valore);
    ftvalt = a_briere .*fT(T,Tmin,Tmax);
    beta  = beta_max*ftvalt*fH(H)*frval;
    dbeta = beta_max*ftvalt*fH(H)*dfrval;
    gammaE = gammaE_max*ftvalt*fH(H)*frval;
    dgammaE = gammaE_max*ftvalt*fH(H)*dfrval;
    muE = muE_min + alphaT_E*abs(T - Topt) + alphaH_E*abs(H - Hopt);
    gammaN = gammaN_max*ftvalt*gW(W);
    muN = muN_min + alphaT_N*abs(T - Topt) + alphaH_N*abs(H - Hopt);
    muN = muN + kW_mort * fW(W);
    muA = muA_min + betaT_A*abs(T - Topt) + betaH_A*abs(H - Hopt);
    % Euler
    E(j+1) = E(j) + ht * ( beta*A(j) - (muE + gammaE)*E(j) );
    N(j+1) = N(j) + ht * ( gammaE*E(j) - (muN + gammaN)*N(j) );
    A(j+1) = A(j) + ht * ( gammaN*N(j) - muA*A(j) );
    sER(j+1) = sER(j) + ht*(dbeta*A(j)-dgammaE*E(j)-(muE+gammaE)*sER(j)+beta*sAR(j)   );
    sNR(j+1) = sNR(j) + ht*(dgammaE*E(j)+gammaE*sER(j)-(muN+gammaN)*sNR(j)  );
    sAR(j+1) = sAR(j) + ht*(gammaN*sNR(j)-muA*sAR(j));
end

%% Plots
sER = sER(1:nday:end);
sNR = sNR(1:nday:end);
sAR = sAR(1:nday:end);

figure
plot(dates,sER,'r','LineWidth',2)
hold on
plot(dates,sNR,'g','LineWidth',2)
hold on
plot(dates,sAR,'b','LineWidth',2)
ylabel('Sensitivity')
title('Sensitivity to $\bar{R}$','Interpreter','latex');
grid on
legend('s_{E,R}','s_{N,R}','s_{A,R}','Location','best')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

% Function Briere
% f_T = (T - Tmin) * sqrt(Tmax - T), if Tmin < T < Tmax; 0 otherwise
function fT_raw = fT(T,Tmin,Tmax)
if T <= Tmin || T >= Tmax
    fT_raw = 0.0;
else
    fT_raw = T*(T - Tmin)*sqrt(max(Tmax - T, 0));
end
end