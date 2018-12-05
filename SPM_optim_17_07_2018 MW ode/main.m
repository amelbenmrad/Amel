clc
clear all
close all
% format long g

%% Experimental values
% importdata Rp_Arash.xlsx
% save Rp_Arash_.mat   ans
load Rp_Arash_.mat
mesures=ans.data;

% Obtenir les valeurs des variables utiles pour ce programme
option=1;   % Rp1
[P,t_exp,Rp_exp] = variables(option,mesures); % This option is run only for 1 second

%%
tfinal=120*60;
CondInit=[P.C1_eq;P.Tb;0;0;0;P.H2i];
global P
[t,X]=ode15s('function_ode_MW',[0 tfinal],CondInit,[]);

%% CI
kp_ref = 38.3305;                            %constante de vitesse de propagation à Tref (m3/mol site/s)
D01= 1.7122e-07;                             %[m²/s]
C1_star=1.0237;
C2_star=   0.7960;

option_opt=0;                                %0: No optimization // 1: optimization
global ic2
if option_opt == 1
    %% Optimisation
    u0=[kp_ref;D01;C1_star;C2_star];
    lb=[20;1e-8;1;0.2];
    ub=[200;1e-6;2;1];
    options_=[];
    u=lsqnonlin(@(x)crit_opt(x,P,t_exp,Rp_exp),u0,lb,ub,options_)
    %C=lsqnonlin(@(x)function(x,param),C0,lb,ub,options)%,[0,0,0,0],ub)%,options)
    %%

    kp_ref = u(1)                            %constante de vitesse de propagation à Tref (m3/mol site/s)
    D01 = u(2)                                %[m²/s]
    C1_star = u(3)
    C2_star = u(4)
end
%% %%%%%%%%%%%%%%%%% BILAN DE MATIERE ET BILAN ENERGETIQUE %%%%%%%%%%%%%%%%%%%%%%

%m peut prendre 0(slab), 1(cylindrical) ou 2(spherical)
m = 2;
%xmesh specifie les points auxquels une solution numérique est demandée pour chaque valeur de tspan (intervalle de concentration)
x = (0: P.r_pol0/100: P.r_pol0);                                        % [m]
%tspan specifie les points auxquels une solution numérique est demandée pour chaque valeur de xmesh (intervalle de temps)
%t = linspace(Start time, number of intervals,End time)
t = 0:0.5:1;                                                            % [s]
%appelle la solution de DiffusionPDEfun1 du composant k pour concentration C1=P au temps tspan(i)
%et au point de maillage xmesh(j) avec CI DiffusionICfun1 et CL DiffusionBCfun1
sol = pdepe(m, @DiffusionPDEfun1, @DiffusionICfun1, @DiffusionBCfun1, x, t, [], P,kp_ref,D01,C1_star,C2_star);

% %%%%%%%%%%%% r_pol %%%%%%%%%%%
res.sol = sol(1:2, :, :);                                               % 3 temps donc sol(1...) c'est t=0 1ere colonne et sol(end...) c'est t=1 derniere colonne
res.t = t(1:2);
P.r_pol = sol(end, :, 4);                                               % [m] 4 pour la solution de r_pol

% %%%%%%%%%%% Rp1_ov %%%%%%%%%%%
C1 = sol(:, :, 1);                                                      % mol/m3(une tranche de la particule)
T = sol(:, :, 3);
kp = kp_ref * exp(-P.Ea./ P.R.*(1 ./ T - 1 / P.T_ref));               % m3/mol/s
kd = P.kd_ref * exp(-P.Ed ./ P.R.*(1 ./ T - 1 / P.T_ref));              % m3/mol/s
C_star = C1_star .* exp(- kd .* t(1:end)') + C2_star;               %(mole site/m3 de cata)
Rp = kp .* C_star .* C1;                                                %(mol/m3/s/dr(m))
dr = [diff(x) , x(end)-x(end-1)];
P.Rp_ov = sum(Rp(3,:).*dr.^3)./sum(dr.^3);                              %(mol/m3(cata)/s) scalaire

Rp1 = Rp * 3600 * P.Mw1 / P.rho_cat;                                    %(g.pol/g.cat/h)
res.Rp1_ov = [sum(Rp1(1,:).*dr.^3)./sum(dr.^3)  ;  sum(Rp1(2,:).*dr.^3)./sum(dr.^3)]; %(g.pol/g.cat/h)

%% Boucle
pas = 20;
for k = 1 : 60*120/pas                                                  % [s]/pas
    ic2= [];    
    ic2.CI2 = sol(end, :, :);                                           % nouvelle condition initiale est la derniere valeur de l'ancien vecteur sol(end..) du precedent    
    ic2.x = (0: P.r_pol/100: P.r_pol);                                  % (m)
    t = (pas*k : pas/2 : pas*(k+1))-9;                                  % [s] k commence à 2 pas a 1 car dans la premiere boucle k-1 = 1
    sol = [];
    sol = pdepe(m, @DiffusionPDEfun1, @DiffusionICfun2, @DiffusionBCfun1, ic2.x, t, [], P,kp_ref,D01,C1_star,C2_star);
    res.sol = [res.sol; sol(1:2, :, :)];                                % 3 temps donc sol(1...) c'est t=0 1ere colonne et sol(end...) c'est t=1 derniere colonne
    res.t = [res.t,t(1:2)];
    r_pol = sol(end, :, 4);                                             % [m] 4 pour la solution de r_pol

    C1 = sol(:, :, 1);                                                  % mol/m3(une tranche de la particule)
    T = sol(:, :, 3);
    kp = kp_ref * exp(-P.Ea./ P.R.*(1 ./ T - 1 / P.T_ref));           % m3/mol/s
    kd = P.kd_ref * exp(-P.Ed ./ P.R.*(1 ./ T - 1 / P.T_ref));          % m3/mol/s
    C_star = C1_star .* exp(- kd .* t(1:end)') + C2_star;           %(mole site/m3 de cata)
    Rp = kp .* C_star .* C1;                                            %(mol/m3/s/dr(m))
    dr = [diff(x) , x(end)-x(end-1)];
    P.Rp_ov = sum(Rp(3,:).*dr.^3)./sum(dr.^3);                          %(mol/m3(cata)/s) scalaire
    
    Rp1 = Rp * 3600 * P.Mw1 / P.rho_cat;                                %(g.pol/g.cat/h)
    res.Rp1_ov = [res.Rp1_ov  ;  sum(Rp1(1,:).*dr.^3)./sum(dr.^3);sum(Rp1(2,:).*dr.^3)./sum(dr.^3)]; %(g.pol/g.cat/h)    

end
Displayplots
