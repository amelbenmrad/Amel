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
option=7;   % Rp1
[P,t_exp,Rp_exp] = variables(option,mesures); % This option is run only for 1 second

%% CI
kp_ref = 70.3;                            %constante de vitesse de propagation à Tref (m3/mol site/s)
D01=2.96e-7;                              %[m²/s]  0.05*

%% Optimisation
u0=[kp_ref;D01];
lb=[50;1e-8];
ub=[150;1e-6];
options_=[];
u=lsqnonlin(@(x)crit_opt(x,P,t_exp,Rp_exp),u0,lb,ub,options_);
%C=lsqnonlin(@(x)function(x,param),C0,lb,ub,options)%,[0,0,0,0],ub)%,options)
Displayplots
