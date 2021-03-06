function [c, f, s] = DiffusionPDEfun1(x, t, u, dudx, P,kp_ref,D01,C1_star,C2_star)
%pdex1pde1 permet de resoudre l'equation partielle diff�rentielle pour le bilan de matiere de l'�thyl�ne
%c(x,t,C,dC/dx)dC/dt = x^-m * d/dx (x^m * f(x,t,C,dC/dx)) + s(x,t,C,dC/dx)

     C1 = u(1);         % Concentration Ethylene [mol/m3 pol]
     C2 = u(2);         % Concentration ICA [mol/m3 pol]
     T = u(3);          % Temperature [K]
     r_pol = u(4);      % Rayon de la particule [m] un scalaire
     nu0 = u(5);        %0th moment of dead chains
     nu1 = u(6);        %1st moment of dead chains
     nu2 = u(7);        %2nd moment of dead chains
     C3 = u(8);         %hydrogen concentration
     
     
    % P.rho_ov � calculer ici (this is variable, depends on [M1] and [M2]!!!!!!!!!!)
    
    kp = kp_ref * exp(-P.Ea / P.R*(1 / T - 1 / P.T_ref));                  %m3/mol/s

    kd = P.kd_ref * exp(-P.Ed / P.R*(1 / T - 1 / P.T_ref));                %m3/mol/s
    
    C_star = C1_star * exp(- kd * t) + C2_star;                            %(mole site/m3 de cata)
   
    mu0 = (P.ki * C_star * C1) / (P.ktH * C3 + kd);                        %0th moment of living chains
     
    mu1 =(P.ki * C_star + kp * mu0) * C1 / (P.ktH * C3 + kd);              %1st moment of living chains 
     
    mu2 = (P.ki * C_star + kp * (2*mu1 + mu0)) * C1 / (P.ktH * C3 + kd);   %2nd moment of living chains
   
    
    Rp = kp * C_star * C1;                                                 %(mol/m3 cata/s)
    
    Rp1 = Rp * 3600 * P.Mw1 / P.rho_cat;                                   %(g.pol/g.cat/h)
    
    phi = r_pol / P.r_cat;                                                 %overall growth factor (-)

    Rv = Rp * ((1-P.epsi) / phi^3);                                        %(mol/m3 cata/s)
   

    %% %%%%%%DIFFUSION%%%%%%

     D_1 = D01*exp(-(P.w1*P.V1_star+P.w2*P.V2_star*P.xi_13/P.xi_23+P.w3*P.V3_star*P.xi_13)/(P.w1*P.K11_div_g*(P.K21_m_Tg+T)+P.w2*P.K12_div_g*(P.K22_m_Tg+T)+P.w3*P.K13_div_g*(P.K23_m_Tg+T)));
%    
%     P.D_2 = P.D02*exp(-(P.w1*P.V1_star*P.xi_23/P.xi_13+P.w2*P.V2_star+P.w3*P.V3_star*P.xi_23)/(P.w1*P.K11_div_g*(P.K21_m_Tg+T)+P.w2*P.K12_div_g*(P.K22_m_Tg+T)+P.w3*P.K13_div_g*(P.K23_m_Tg+T)));

    %% %%%%% PDE COEFFICIENT %%%%%
%       C1                  C2                  T                       r_pol                                                   nu0                      nu1                      nu2                 C3  
    c = [1;                 1;                  P.rho_ov * P.Cp_ov;     1;                                                      1;                       1;                       1;                  1];
    f = [D_1 * dudx(1);     P.D_2 * dudx(2);    P.kc_ov * dudx(3);      0;                                                      0;                       0;                       0;                  0];
    s = [-Rv;               0;                  -P.delta_Hp * Rv;       ((1-P.epsi) / phi^3)*P.Rp_ov*r_pol*P.Mw1/3/P.rho_ov;    (P.ktH*C3+kd)*mu0;       (P.ktH*C3+kd)*mu1;       (P.ktH*C3+kd)*mu2;   P.ktH*C3*mu0];

end