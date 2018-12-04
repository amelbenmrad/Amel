
function [pl, ql, pr, qr] = DiffusionBCfun1(xl, ul, xr, ur, t, P,kp_ref,D01)
%pdex1bc1 permet d'implémenter des conditions limite à l'equation partielle différentielle de l'éthylène
%p(x,t,C) + q(x,t) * f(x,t,C,dC/dx) = 0
    
    %% %%%%%%HEAT TRANSFER COEFFICIENT%%%%%%
    T = ur(3);
    r_pol = ur(4);
    
    rho_g = P.P * P.Mw1 / (P.R * T);                                    %densité de phase gaz (kg/m3)
          
    %%%Nicollela correllation%%%
%      Pr = (P.mu_g * P.Cp_g) / P.kc_g;
%      h = (P.kc_g / (2 * r_pol)) * (2 + 0.265 * ((2 * r_pol * P.epsi_s * rho_g^3)/(P.mu_g^3)^0.241) * Pr^(1/3));    %(J/m2/s/K)

    %%%Ranz-Marshall correlation%%%
       %P.v = 1 - sin(t/2);
    Re = (rho_g * P.v * (2 * r_pol)) / P.mu_g;
    Pr = (P.mu_g * P.Cp_g) / P.kc_g;
    h = (P.kc_g / (2 * r_pol)) * (2 + 0.6 * Re^(1/2) * Pr^(1/3));                                                    %(J/m2/s/K)

    %% %%%%%%DIFFUSION%%%%%%
    
     D_1=D01*exp(-(P.w1*P.V1_star+P.w2*P.V2_star*P.xi_13/P.xi_23+P.w3*P.V3_star*P.xi_13)/(P.w1*P.K11_div_g*(P.K21_m_Tg+T)+P.w2*P.K12_div_g*(P.K22_m_Tg+T)+P.w3*P.K13_div_g*(P.K23_m_Tg+T)));
    
    P.D_2=P.D02*exp(-(P.w1*P.V1_star*P.xi_23/P.xi_13+P.w2*P.V2_star+P.w3*P.V3_star*P.xi_23)/(P.w1*P.K11_div_g*(P.K21_m_Tg+T)+P.w2*P.K12_div_g*(P.K22_m_Tg+T)+P.w3*P.K13_div_g*(P.K23_m_Tg+T)));

    %% %%%%% PDE COEFFICIENT %%%%%
    %    [M1;                 M2,                   T                   r_pol]
    pl = [0;                  0;                    0   ;               0];
    ql = [1 / D_1;          1 / P.D_2;            1 / P.kc_ov;        1];
    pr = [ur(1) - P.C1_eq;    ur(2) - P.C2_eq;      h * (ur(3) - P.Tb); 0];
    qr = [0;                  0;                    1;                  1];
    
end