clc
clear all
close all

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

%% Plots
figure
hold all
plot(t/60, X(:,1), 'LineWidth', 2) 
xlabel('Time (min)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel('(mol/m^3)', 'fontsize', 12, 'fontweight' , 'b', 'fontname', 'arial')
% axis([0 inf 20 inf])
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')
title('Monomer concentration')

figure
hold all
plot(t/60, X(:,6), 'LineWidth', 2) 
xlabel('Time (min)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel(' (mol/m^3)', 'fontsize', 12, 'fontweight' , 'b', 'fontname', 'arial')
% axis([0 inf 20 inf])
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')
title('Hydrogen concentration')

figure
hold all
plot(t/60, X(:,3:5), 'LineWidth', 2) 
xlabel('Time (min)', 'fontsize', 12, 'fontweight', 'b', 'fontname', 'arial')
ylabel(' (mol/m^3)', 'fontsize', 12, 'fontweight' , 'b', 'fontname', 'arial')
% axis([0 inf 20 inf])
set(gca, 'FontSize', 12, 'fontweight', 'b', 'fontname', 'arial')
title('Moments of dead polymer')