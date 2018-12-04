function [J] = crit_opt(u,P,t_exp,Rp_exp)

global ic2
kp_ref = u(1)                            %constante de vitesse de propagation à Tref (m3/mol site/s)
D01= u(2)                                %diffusion [m²/s]

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
sol = pdepe(m, @DiffusionPDEfun1, @DiffusionICfun1, @DiffusionBCfun1, x, t, [], P,kp_ref,D01);
 
% %%%%%%%%%%%% r_pol %%%%%%%%%%%
res.sol = sol(1:2, :, :);                                               % 3 temps donc sol(1...) c'est t=0 1ere colonne et sol(end...) c'est t=1 derniere colonne
res.t = t(1:2);
P.r_pol = sol(end, :, 4);                                               % [m] 4 pour la solution de r_pol

% %%%%%%%%%%% Rp1_ov %%%%%%%%%%%
C1 = sol(:, :, 1);                                                      % mol/m3(une tranche de la particule)
T = sol(:, :, 3);
kp = kp_ref * exp(-P.Ea./ P.R.*(1 ./ T - 1 / P.T_ref));                 % m3/mol/s
kd = P.kd_ref * exp(-P.Ed ./ P.R.*(1 ./ T - 1 / P.T_ref));              % m3/mol/s
C_star = P.C1_star .* exp(- kd .* t(1:end)') + P.C2_star;               %(mole site/m3 de cata)
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
    sol = pdepe(m, @DiffusionPDEfun1, @DiffusionICfun2, @DiffusionBCfun1, ic2.x, t, [], P,kp_ref,D01);
    res.sol = [res.sol; sol(1:2, :, :)];                                % 3 temps donc sol(1...) c'est t=0 1ere colonne et sol(end...) c'est t=1 derniere colonne
    res.t = [res.t,t(1:2)];
    r_pol = sol(end, :, 4);                                             % [m] 4 pour la solution de r_pol

    C1 = sol(:, :, 1);                                                  % mol/m3(une tranche de la particule)
    T = sol(:, :, 3);
    kp = kp_ref * exp(-P.Ea./ P.R.*(1 ./ T - 1 / P.T_ref));             % m3/mol/s
    kd = P.kd_ref * exp(-P.Ed ./ P.R.*(1 ./ T - 1 / P.T_ref));          % m3/mol/s
    C_star = P.C1_star .* exp(- kd .* t(1:end)') + P.C2_star;           %(mole site/m3 de cata)
    Rp = kp .* C_star .* C1;                                            %(mol/m3/s/dr(m))
    dr = [diff(x) , x(end)-x(end-1)];
    P.Rp_ov = sum(Rp(3,:).*dr.^3)./sum(dr.^3);                          %(mol/m3(cata)/s) scalaire
    
    Rp1 = Rp * 3600 * P.Mw1 / P.rho_cat;                                %(g.pol/g.cat/h)
    res.Rp1_ov = [res.Rp1_ov  ;  sum(Rp1(1,:).*dr.^3)./sum(dr.^3);sum(Rp1(2,:).*dr.^3)./sum(dr.^3)]; %(g.pol/g.cat/h)    

end

%% Critère pour l'optimisation
for k1 = 1:length(t_exp)
    %ind1(k1) = min(find(abs(t_exp(k1)*60 - res.t) < 2.5));
    ind1(k1) = min(find(abs(t_exp(k1)*60 - res.t) < 5));
end

J = (abs(res.Rp1_ov(ind1) - Rp_exp));
sum(abs(J))

%% Affichage optimisation
figure(1)
clf
hold all
plot(res.t/60, res.Rp1_ov, 'LineWidth', 2)   %sum(A,2) = vecteur colonne contenant la somme des lignes du vecteur A
hold on
plot(t_exp, Rp_exp, 'm*')   %sum(A,2) = vecteur colonne contenant la somme des lignes du vecteur A
legend('Model','Experimental')
xlabel('Time (min)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel('Reaction rate (g.pol/g.cat/h)', 'fontsize', 12, 'fontweight' , 'b', 'fontname', 'arial')
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')
pause(0.0001)