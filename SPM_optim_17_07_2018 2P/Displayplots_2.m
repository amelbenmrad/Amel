% sol = t,r,etat(C1,C2, T, r_pol)
C1 = res.sol(:, :, 1);
C2 = res.sol(:, :, 2);
T = res.sol(:, :, 3);
r_pol = res.sol(:, end, 4); % [m]

for k=1:length(r_pol)
    x(k,:) = (0: r_pol(k)/100: r_pol(k));                                        % [m]
end

kp = P.kp_ref * exp(-P.Ea./ P.R.*(1 ./ T - 1 / P.T_ref));        % m3/mol/s

kd = P.kd_ref * exp(-P.Ed ./ P.R.*(1 ./ T - 1 / P.T_ref));       % m3/mol/s

C_star = P.C1_star .* exp(- kd .* res.t') + P.C2_star;           %(mole site/m3 de cata)

Rp = kp .* C_star .* C1;                                         %(mol/m3/s)

Rp1 = Rp * 3600 * P.Mw1 / P.rho_cat;                             %(g.pol/g.cat/h)

phi = r_pol ./ P.r_cat;                                          %overall growth factor (-)

Rv = Rp .* ((1-P.epsi) ./ phi.^3);                               %(mol/m3/s)
  
%% %%%CATALYST CONCENTRATION%%%

figure
hold all
plot(res.t/60, C_star, 'LineWidth', 2)   %sum(A,2) = vecteur colonne contenant la somme des lignes du vecteur A
xlabel('Time (min)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel('Catalyst concentration (mol/m^3)', 'fontsize', 12, 'fontweight' , 'b', 'fontname', 'arial')
% axis([0 inf 20 inf])
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')
title('Catalyst concentration')

%% %%%REACTION RATE%%%
figure
hold all
plot(res.t/60, Rp1, 'LineWidth', 2)   %sum(A,2) = vecteur colonne contenant la somme des lignes du vecteur A
xlabel('Time (min)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel('Reaction rate (g.pol/g.cat/h)', 'fontsize', 12, 'fontweight' , 'b', 'fontname', 'arial')
% axis([0 inf 0 inf])
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')

figure
hold all
plot(res.t/60, res.Rp1_ov, 'LineWidth', 2)   %sum(A,2) = vecteur colonne contenant la somme des lignes du vecteur A
hold on
plot(t_exp, Rp_exp, 'm*')   %sum(A,2) = vecteur colonne contenant la somme des lignes du vecteur A
legend('Model','Experimental')
xlabel('Time (min)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel('Reaction rate (g.pol/g.cat/h)', 'fontsize', 12, 'fontweight' , 'b', 'fontname', 'arial')
% Saxis([0 inf 0 1700])
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')
 
% figure
% hold all
% plot(res.t/60, Rp1(:, 90), 'LineWidth', 2)   %sum(A,2) = vecteur colonne contenant la somme des lignes du vecteur A
% xlabel('Time (min)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
% ylabel('Reaction rate (g.pol/g.cat/h)', 'fontsize', 12, 'fontweight' , 'b', 'fontname', 'arial')
% % axis([0 inf 0 inf])
% set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')
%% Radius
figure
hold all
[i1,i2]=min(abs(res.t-300));plot(x(i2,:)/x(i2,end) ,Rp1(i2, :) , 'LineWidth', 2)
[i1,i2]=min(abs(res.t-600));plot(x(i2,:)/x(i2,end),Rp1(i2, :) , 'k--', 'LineWidth', 2)
[i1,i2]=min(abs(res.t-1200));plot(x(i2,:)/x(i2,end) ,Rp1(i2, :) , 'm:', 'LineWidth', 2)
[i1,i2]=min(abs(res.t-2400));plot(x(i2,:)/x(i2,end),Rp1(i2, :), 'g:', 'LineWidth', 2)
legend('t = 5 min', 't = 10 min', 't = 20 min', 't = 40 min')
xlabel('Normalized particle radius (-)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel('Reaction rate (g.pol/g.cat/h)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
% axis([0 inf 80 90])
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')

%% %%%PARTICLE RADIUS%%%

figure
hold all 
plot(res.t/60, r_pol * 10^6, 'LineWidth', 2)   %sum(A,2) = vecteur colonne contenant la somme des lignes du vecteur A
%res.t/60 pour avoir 1 min en temps dans l'axe des x
xlabel('Time (min)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel('Particle radius (µm)', 'fontsize', 12, 'fontweight' , 'b', 'fontname', 'arial')
% axis([0 inf 20 inf])
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')

v_pol = 4 / 3 * pi() * (r_pol.^3 - P.r_pol0^3*0);
m_pol=v_pol * 869.8 * 1000; % (g)

figure
hold all
plot(res.t/60, m_pol, 'LineWidth', 2)   %sum(A,2) = vecteur colonne contenant la somme des lignes du vecteur A
xlabel('Time (min)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel('Total mass of polymer (g)', 'fontsize', 12, 'fontweight' , 'b', 'fontname', 'arial')
% axis([0 inf 20 inf])
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')

%% %%%CONCENTRATION ETHYLENE%%%%%

% % 3D surface plot
% figure
% surf(x*1e6,t,C1,'edgecolor','none');
% xlabel('Particle radius (µm)','fontsize',12,'fontweight','b','fontname','arial')
% ylabel('Time (s)','fontsize',12,'fontweight','b','fontname','arial')
% zlabel('Ethylene concentration (mol m^{-3})','fontsize',12,'fontweight','b','fontname','arial')
% axis([0 inf 0 inf 0 inf])
% set(gcf(), 'Renderer', 'painters')
% set(gca,'FontSize',12,'fontweight','b','fontname','arial')

figure
hold all
[i1,i2]=min(abs(res.t-300));plot(x(i2,:)/x(i2,end) , C1(i2, :), 'LineWidth', 2)        
[i1,i2]=min(abs(res.t-600));plot(x(i2,:)/x(i2,end) , C1(i2, :), 'k--', 'LineWidth', 2)
[i1,i2]=min(abs(res.t-1200));plot(x(i2,:)/x(i2,end) , C1(i2, :), 'm:', 'LineWidth', 2)
[i1,i2]=min(abs(res.t-2400));plot(x(i2,:)/x(i2,end) , C1(i2, :), 'g', 'LineWidth', 2)
legend('t = 5 min', 't = 10 min', 't = 20 min', 't = 40 min')
xlabel('Normalized particle radius (-)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel('Ethylene concentration (mol m^{-3})', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
axis([0 inf 0 200])
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')

figure
hold all
plot(res.t/60, C1(:, end), 'LineWidth', 2)
plot(res.t/60, C1(:, round(length(x(1,:)) / 2)), 'k--', 'LineWidth', 2)            % round = arrondi au plus haut
plot(res.t/60, C1(:, 1), 'm:', 'LineWidth', 2)
legend('R', 'R/2', 'Center')
xlabel('Time (min)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel('Ethylene concentration (mol m^{-3})', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
% axis([0 inf 20 inf])
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')

%% %%%CONCENTRATION ICA%%%%%

% 3D surface plot
% figure
% surf(x*1e6,t,C2,'edgecolor','none');
% xlabel('Particle radius (µm)','fontsize',12,'fontweight','b','fontname','arial')
% ylabel('Time (s)','fontsize',12,'fontweight','b','fontname','arial')
% zlabel('ICA concentration (mol m^{-3})','fontsize',12,'fontweight','b','fontname','arial')
% axis([0 inf 0 inf 0 inf])
% set(gcf(), 'Renderer', 'painters')
% set(gca,'FontSize',12,'fontweight','b','fontname','arial')

figure
hold all
[i1,i2]=min(abs(res.t-300));plot(x(i2,:)/x(i2,end) , C2(i2, :), 'LineWidth', 2)
[i1,i2]=min(abs(res.t-600));plot(x(i2,:)/x(i2,end) , C2(i2, :), 'k--', 'LineWidth', 2)
[i1,i2]=min(abs(res.t-1200));plot(x(i2,:)/x(i2,end) , C2(i2, :), 'm:', 'LineWidth', 2)
[i1,i2]=min(abs(res.t-2400));plot(x(i2,:)/x(i2,end) , C2(i2, :) ,'g', 'LineWidth', 2)
legend('t = 5 min', 't = 10 min', 't = 20 min', 't = 40 min')
xlabel('Normalized particle radius (-)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel('ICA concentration (mol m^{-3})', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
%axis([0 inf 0 inf])
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')

figure
hold all
plot(res.t/60, C2(:, end), 'LineWidth', 2)   %Linewidth = largeur de la ligne
plot(res.t/60, C2(:, round(length(x(1,:)) / 2)), 'k--', 'LineWidth', 2) % round = arrondi au plus haut
plot(res.t/60, C2(:, 1), 'm:','LineWidth', 2)
legend('R', 'R/2', 'Center')
xlabel('Time (min)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel('ICA concentration (mol m^{-3})', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
% axis([0 inf 20 inf])
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')

%% %%%TEMPERATURE%%%%%

% 3D surface plot
% figure
% surf(x*1e6,t,T,'edgecolor','none');
% xlabel('Particle radius x (µm)','fontsize',12,'fontweight','b','fontname','arial')
% ylabel('Time (s)','fontsize',12,'fontweight','b','fontname','arial')
% zlabel('Temperature (K)','fontsize',12,'fontweight','b','fontname','arial')
% axis([0 inf 0 inf 0 inf])
% set(gcf(), 'Renderer', 'painters')
% set(gca,'FontSize',12,'fontweight','b','fontname','arial')

figure
hold all
[i1,i2]=min(abs(res.t-300));plot(x(i2,:)/x(i2,end) ,T(i2, :), 'LineWidth', 2)
[i1,i2]=min(abs(res.t-600));plot(x(i2,:)/x(i2,end) ,T(i2, :), 'k--', 'LineWidth', 2)
[i1,i2]=min(abs(res.t-1200));plot(x(i2,:)/x(i2,end) ,T(i2, :) , 'm:', 'LineWidth', 2)
[i1,i2]=min(abs(res.t-2400));plot(x(i2,:)/x(i2,end) ,T(i2, :) , 'g:', 'LineWidth', 2)
%plot(x * 1e6,T(res.t == 500, :) , 'b:', 'LineWidth', 2)
legend('t = 5 min', 't = 10 min', 't = 20 min', 't = 40 min')
xlabel('Normalized particle radius (-)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel('Temperature (°C)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
% axis([0 inf 80 90])
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')

figure
hold all
plot(res.t/60, T(:, end) - 273, 'LineWidth' ,2)
plot(res.t/60, T(:, round(length(x(1,:))/2)) - 273, 'k--', 'LineWidth', 2)
plot(res.t/60, T(:, 1) - 273, 'm:', 'LineWidth', 2)
legend('Surface', 'Surface/2', 'Center')
xlabel('Time (min)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel('Temperature (°C)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
% axis([0 inf 20 inf])
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')

figure
hold all
plot(res.t/60, T(:, end) -  P.Tb, 'LineWidth' ,2)
plot(res.t/60, T(:, round(length(x(1,:))/2)) -  P.Tb, 'k--', 'LineWidth', 2)
plot(res.t/60, T(:, 1) -  P.Tb, 'm:', 'LineWidth', 2)
legend('Surface', 'Surface/2', 'Center')
xlabel('Time (min)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel('T_s - T_b (°C)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
axis([0 inf 0 20])
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')