function dXdt = function_ode_MW(t,u,P)
global P
     C1 = u(1);         % Concentration Ethylene [mol/m3 pol]
%    C2 = u(2);         % Concentration ICA [mol/m3 pol]
     T = u(2);          % Temperature [K]
%    r_pol = u(4);      % Rayon de la particule [m] un scalaire
     nu0 = u(3);        %0th moment of dead chains
     nu1 = u(4);        %1st moment of dead chains
     nu2 = u(5);        %2nd moment of dead chains
     C3 = u(6);         %hydrogen concentration
     
   
    kp = P.kp_ref * exp(-P.Ea / P.R*(1 / T - 1 / P.T_ref));                 % m3/mol/s

    kd = P.kd_ref * exp(-P.Ed / P.R*(1 / T - 1 / P.T_ref));                 % m3/mol/s
    
    C_star = P.C1_star * exp(- kd * t) + P.C2_star;                         %(mole site/m3 de cata)
   
    mu0 = (P.ki * C_star * C1) / (P.ktH * C3 + kd) ;                        %0th moment of living chains
     
    mu1 =(P.ki * C_star + kp * mu0) * C1 / (P.ktH * C3 + kd);               %1st moment of living chains 
     
    mu2 = (P.ki * C_star + kp * (2*mu1 + mu0)) * C1 / (P.ktH * C3 + kd);    %2nd moment of living chains
   
    Rp = kp * C_star * C1;                                                  %(mol/m3 cata/s)
    
%     Rp1 = Rp * 3600 * P.Mw1 / P.rho_cat;                                  %(g.pol/g.cat/h)

    phi = 0.9;% r_pol / P.r_cat;                                            %overall growth factor (-)

    Rv = Rp * ((1-P.epsi) / phi^3);                                         %(mol/m3 cata/s)

%% differential equations
%========================
F=1000;             %monomer flow rate (mol/m3 cata/s)
dC1dt=F-Rv;
dTdt=0;             %-P.delta_Hp * Rv ;
dnu0dt=(P.ktH*C3+kd)*mu0;
dnu1dt=(P.ktH*C3+kd)*mu1;
dnu2dt=(P.ktH*C3+kd)*mu2;
Fc3=1;              %hydrogen flow rate (mol/s) ??
dC3dt=Fc3-P.ktH*C3*mu0;

dXdt=[dC1dt;dTdt;dnu0dt;dnu1dt;dnu2dt;dC3dt];
end