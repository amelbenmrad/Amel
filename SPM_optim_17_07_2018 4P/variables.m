function [P,t_exp,Rp_exp] = variables(option,mesures)

    % Declare et initialise toutes les variables utilisées dans le programme principal "main"
    
    %% %%POLYMERIZATION KINETICS PARAMETERS%%%%
    %P.kp_ref = 70.3;                            %constante de vitesse de propagation à Tref (m3/mol site/s)
    P.kd_ref = 1.56e-4;                         %constante de desactivation (s-1)
%     P.C1_star = 1.37;                          %constante déterminant concentration aux SA au début de réaction (mole site/m3 de cata)
%     P.C2_star = 0.45;                           %constante déterminant concentration aux SA à la fin de réaction (mole site/m3 de cata)
    P.Ea = 42000;                               %Energie d'activation (J/mol)
    P.Ed = 42000;                               %Energie de desactivation (J/mol)
    P.T_ref = 313.15;                           %température de réference (K)
    P.delta_Hp = -107500;                       %enthalpie de polymérisation de l'éthylène (J/mol)
    
    %% %%HYDRODYNAMICS OF REACTOR%%%%
    P.v = 0.02;                                 %vitesse relative gaz-particule (m/s) (Arash 0.02, Kosek 0.2)
    % pour étudier des vitesses différentes: v = 1 ? sin(t/2)
    
    %% %%CATALYST PHYSICAL PROPERTIES%%%%
    P.r_cat = 25e-6;                            %rayon de particule de catalyseur initiale (m)
    P.rho_cat = 2300;                           %densité du catalyseur (kg/m3)
    P.R = 8.314;                                %constante des gaz parfaits (m3.Pa/mol/K)
    
    %% %%OVERALL PARTICLE PROPERTIES%%%%
    P.r_pol0 = 50e-6;                           %rayon de la particule de polymere seulement (sans le volume des pores) (m) (dp=100 micrometre)
    P.epsi = 0.05;                              %porosité de la particule (-)
    P.kc_ov = 0.2;                              %conductivité thermique de la particule (J/m.s.K) (Kiparissides 0.108, Ray, et Arash =0.2)
    P.Cp_ov = 1840;                             %capacité calorifique de la particule (J/kg.K)
    P.epsi_s = 0.4;                            %viscous dissipation rate (m2/s3) (A VERIFIER !!!!!!!)
    
    %% %%GAZ PHASE PROPERTIES%%%%
    P.Tb = 353.15;                              %température bulk phase (K)
    P.Mw1 = 0.028;                              %molecular weight of ethylene (kg/mol)
    P.Mw2=0.06818;                              % molecular weight of n-hexane (kg/mol)
    
    P.Rp_ov=0;
    
    
    %% Diffusion parameters
    %P.D01=2.96e-7;              % [m²/s]  0.05*
    P.D02=3.5e-8;               % [m²/s]
    
    P.V1_star=1.341;            % [cm3/g]
    P.V2_star=1.133;            % [cm3/g]
    P.V3_star=1.006;            % [cm3/g]
    
    P.K11_div_g=1.97e-3;        % [cm3/g.K]
    P.K12_div_g=1.96e-3;        % [cm3/g.K]
    P.K13_div_g=1.02e-3;        % [cm3/g.K]
    
    P.K21_m_Tg=42.38;           % [K]
    P.K22_m_Tg=-41.08;          % [K]
    P.K23_m_Tg=-228.7;          % [K]
    
    P.xi_13=0.4548;             % [-]
    P.xi_23=0.9184;             % [-]
   
    %% %%OPTIONS%%%%
    if option==1                                %P_eth=7bar P_ICA=0.3bar
        %%OVERALL PARTICLE PROPERTIES%%%%
        P.C1_eq = 86; %69.39;                    %concentration à l'éq de l'éthylène dans particule (mol/m3) comes from Sanchez Lacombe
        P.C2_eq = 122.6; %123.85;                %concentration du n-hexane dans particule (mol/m3) comes from Sanchez Lacombe
        P.rho_ov = 869.8;                        %densité de particule (kg/m3) (this is variable, depnds on [M1] and [M2]!!!!!!!!!!)
        %P.D_1 = 3.1e-9;                         %diffusivité totale de l'éthylène (m2/s)
        P.D_2 = 5.2e-10;                        %diffusivité totale du n-hexane (m2/s)
        P.w1 = 0.00583;
        P.w2 = 0.0375;
        P.w3=1-P.w1-P.w2;
        
        %%GAZ PHASE PROPERTIES%%%%
        P.P = 7e5;                              %pressure of ethylene(Pa)
        P.mu_g = 1.32e-5;                       %viscosité de la phase gaz (kg/m.s)
        P.kc_g = 2.61e-2;                       %conductivité thermique de phase gaz (J/m.s.K)
        P.Cp_g =1606;                           %capacité calorifique de la phase gaz (J/kg.K)
        %P.rho_g = 9.2;                         %densité de phase gaz (kg/m3)
          
    elseif option==2                            %P_eth=7bar P_ICA=0.6bar
        %%OVERALL PARTICLE PROPERTIES%%%%
        P.C1_eq = 93.1; %73.35;                  %concentration à l'éq de l'éthylène dans particule (mol/m3) comes from Sanshez Lacomb
        P.C2_eq = 264.9; %262.69;                %concentration du n-hexane dans particule (mol/m3)
        P.rho_ov = 864.4;                        %densité de particule (kg/m3) (this is variable, depnds on [M1] and [M2]!!!!!!!!!!)
        P.D_1 = 4.4e-9;                         %diffusivité totale de l'éthylène (m2/s)
        P.D_2 = 5e-10;                          %diffusivité totale du n-hexane (m2/s)
        P.w1 = 0.00739;
        P.w2 = 0.08126;
        P.w3=1-P.w1-P.w2;
        
        %%GAZ PHASE PROPERTIES%%%%
        P.P = 7e5;                              %pressure of ethylene(Pa)
        P.mu_g = 1.29e-5;                       %viscosité de la phase gaz (kg/m.s)
        P.kc_g = 2.56e-2;                       %conductivité thermique de phase gaz (J/m.s.K)
        P.Cp_g =1642;                           %capacité calorifique de la phase gaz (J/kg.K)
        %P.rho_g = 10.9;                        %densité de phase gaz (kg/m3)
        
    elseif option==3                            %P_eth=7bar P_ICA=0.8bar NON
        %%OVERALL PARTICLE PROPERTIES%%%%
        P.C1_eq = 99.1; %76.64;                  %concentration à l'éq de l'éthylène dans particule (mol/m3) comes from Sanshez Lacomb
        P.C2_eq = 381.5; %373.99;                %concentration du n-hexane dans particule (mol/m3)
        P.rho_ov = 859.9;                        %densité de particule (kg/m3) (this is variable, depnds on [M1] and [M2]!!!!!!!!!!)
        P.D_1 = 7.8e-9;                         %diffusivité totale de l'éthylène (m2/s)
        P.D_2 = 4.9e-10;                        %diffusivité totale du n-hexane (m2/s)
        P.w1 = 0.007876;
        P.w2 = 0.118;
        P.w3=1-P.w1-P.w2;
        
        %%GAZ PHASE PROPERTIES%%%%
        P.P = 7e5;                              %pressure of ethylene(Pa)
        P.mu_g = 1.28e-5;                       %viscosité de la phase gaz (kg/m.s)
        P.kc_g = 2.52e-2;                       %conductivité thermique de phase gaz (J/m.s.K)
        P.Cp_g =1664;                           %capacité calorifique de la phase gaz (J/kg.K)
        %P.rho_g = 10.9;                        %densité de phase gaz (kg/m3)
        
    elseif option==4                            %%P_eth=12bar P_ICA=0.3bar
        %%OVERALL PARTICLE PROPERTIES%%%%
        P.C1_eq = 146.6;                        %concentration à l'éq de l'éthylène dans particule (mol/m3) comes from Sanshez Lacomb
        P.C2_eq = 116.8;                        %concentration du n-hexane dans particule (mol/m3)
        P.rho_ov = 868.7;                       %densité de particule (kg/m3) (this is variable, depnds on [M1] and [M2]!!!!!!!!!!)
        P.D_1 = 1.7e-9;                         %diffusivité totale de l'éthylène (m2/s)
        P.D_2 = 3.2e-10;                        %diffusivité totale du n-hexane (m2/s)
        %P.w1 = 0.01163;
        %P.w2 = 0.03555;
        %P.w3=1-P.w1-P.w2;

        %%GAZ PHASE PROPERTIES%%%%
        P.P = 12e5;                             %pressure of ethylene(Pa)
        P.mu_g = 1.30e-5;                       %viscosité de la phase gaz (kg/m.s)
        P.kc_g = 2.67e-2;                       %conductivité thermique de phase gaz (J/m.s.K)
        P.Cp_g =1697;                           %capacité calorifique de la phase gaz (J/kg.K)
        %P.rho_g = 14.5;                        %densité de phase gaz (kg/m3)
        
    elseif option==5                            %%P_eth=12bar P_ICA=0.6bar NON
        %%OVERALL PARTICLE PROPERTIES%%%%
        P.C1_eq = 158.5;                        %concentration à l'éq de l'éthylène dans particule (mol/m3) comes from Sanshez Lacomb
        P.C2_eq = 259.2;                        %concentration du n-hexane dans particule (mol/m3)
        P.rho_ov = 863.1;                       %densité de particule (kg/m3) (this is variable, depnds on [M1] and [M2]!!!!!!!!!!)
        %P.D_1 = 2e-9;                         %diffusivité totale de l'éthylène (m2/s)
        %P.D_2 = 3.1e-10;                        %diffusivité totale du n-hexane (m2/s)
        P.w1 = 0.012556;
        P.w2 = 0.07923;
        P.w3=1-P.w1-P.w2;
        
        %%GAZ PHASE PROPERTIES%%%%
        P.P = 12e5;                             %pressure of ethylene(Pa)
        P.mu_g = 1.28e-5;                       %viscosité de la phase gaz (kg/m.s)
        P.kc_g = 2.63e-2;                       %conductivité thermique de phase gaz (J/m.s.K)
        P.Cp_g =1721;                           %capacité calorifique de la phase gaz (J/kg.K)
        %P.rho_g = 14.5;                        %densité de phase gaz (kg/m3)
        
    elseif option==6                                        %P_eth=12bar P_ICA=0.8bar NON
        %%OVERALL PARTICLE PROPERTIES%%%%
        P.C1_eq = 168.8;                        %concentration à l'éq de l'éthylène dans particule (mol/m3) comes from Sanshez Lacomb
        P.C2_eq = 379.8;                        %concentration du n-hexane dans particule (mol/m3)
        P.rho_ov = 858.3;                       %densité de particule (kg/m3) (this is variable, depnds on [M1] and [M2]!!!!!!!!!!)
        %P.D_1 = 3.3e-9;                         %diffusivité totale de l'éthylène (m2/s)
        %P.D_2 = 3.1e-10;                        %diffusivité totale du n-hexane (m2/s)
        P.w1 = 0.01339;
        P.w2 = 0.1171;
        P.w3=1-P.w1-P.w2;
        
        %%GAZ PHASE PROPERTIES%%%%
        P.P = 12e5;                             %pressure of ethylene(Pa)
        P.mu_g = 1.27e-5;                       %viscosité de la phase gaz (kg/m.s)
        P.kc_g = 2.61e-2;                       %conductivité thermique de phase gaz (J/m.s.K)
        P.Cp_g =1735;                           %capacité calorifique de la phase gaz (J/kg.K)
        %P.rho_g = 16.2;                        %densité de phase gaz (kg/m3)
        
    elseif option==7                             %P_eth=7bar without ICA (Rp1)
        %%OVERALL PARTICLE PROPERTIES%%%%
        P.C1_eq = 80.2;                         %concentration à l'éq de l'éthylène dans particule (mol/m3) comes from Sanchez Lacombe
        P.C2_eq = 0;                            %concentration du n-hexane dans particule (mol/m3) comes from Sanchez Lacombe
        P.rho_ov = 874.5;                       %densité de particule (kg/m3) (this is variable, depnds on [M1] and [M2]!!!!!!!!!!)
%         P.D_1 = 2.6e-9;                         %diffusivité totale de l'éthylène (m2/s)
        P.D_2 = 0;                              %diffusivité totale du n-hexane (m2/s)
       P.w1 = 0.01339; % FAUX
        P.w2 = 0.1171; % FAUX
        P.w3=1-P.w1-P.w2;
               
        %%GAZ PHASE PROPERTIES%%%%
        P.P = 7e5;                              %pressure of ethylene(Pa)
        P.mu_g = 1.36e-5;                       %viscosité de la phase gaz (kg/m.s)
        P.kc_g = 2.68e-2;                       %conductivité thermique de phase gaz (J/m.s.K)
        P.Cp_g =1563;                           %capacité calorifique de la phase gaz (J/kg.K)
        %P.rho_g = 8.3;                         %densité de phase gaz (kg/m3)

    else                                        %P_eth=12bar without ICA (Rp5)
        %%OVERALL PARTICLE PROPERTIES%%%%
        P.C1_eq =137.3;                         %concentration à l'éq de l'éthylène dans particule (mol/m3) comes from Sanchez Lacombe
        P.C2_eq =0;                             %concentration du n-hexane dans particule (mol/m3) comes from Sanchez Lacombe
        P.rho_ov = 873.2;                       %densité de particule (kg/m3) (this is variable, depnds on [M1] and [M2]!!!!!!!!!!)
        P.D_1 = 1.4e-9;                         %diffusivité totale de l'éthylène (m2/s)
        P.D_2 = 0;                              %diffusivité totale du n-hexane (m2/s)
%         P.w1 =
%         P.w2 = 
%         P.w3=1-P.w1-P.w2;
        
        %%GAZ PHASE PROPERTIES%%%%
        P.P = 12e5;                              %pressure of ethylene(Pa)
        P.mu_g = 1.31e-5;                       %viscosité de la phase gaz (kg/m.s)
        P.kc_g = 2.71e-2;                       %conductivité thermique de phase gaz (J/m.s.K)
        P.Cp_g =1672;                           %capacité calorifique de la phase gaz (J/kg.K)
        %P.rho_g = 13.4;                         %densité de phase gaz (kg/m3)
    end
    if option==1   
    t_exp=mesures.Rp2(1:end,1); % [min]
    Rp_exp=mesures.Rp2(1:end,2); % g polym / g cata / h
elseif option==2
    t_exp=mesures.Rp3(1:end,1); % [min]
    Rp_exp=mesures.Rp3(1:end,2); % g polym / g cata / h
elseif option==3
    t_exp=mesures.Rp4(1:end,1); % [min]
    Rp_exp=mesures.Rp4(1:end,2); % g polym / g cata / h
elseif option==4
    t_exp=mesures.Rp6(1:end,1); % [min]
    Rp_exp=mesures.Rp6(1:end,2); % g polym / g cata / h
elseif option==5
    t_exp=mesures.Rp7(1:end,1); % [min]
    Rp_exp=mesures.Rp7(1:end,2); % g polym / g cata / h
elseif option==6
    t_exp=mesures.Rp8(1:end,1); % [min]
    Rp_exp=mesures.Rp8(1:end,2); % g polym / g cata / h
elseif option==7
    t_exp=mesures.Rp1(1:end,1); % [min]
    Rp_exp=mesures.Rp1(1:end,2); % g polym / g cata / helseif option==6
elseif option==8
    t_exp=mesures.Rp5(1:end,1); % [min]
    Rp_exp=mesures.Rp5(1:end,2); % g polym / g cata / h
end    
end