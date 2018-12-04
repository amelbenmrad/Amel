function u0 = DiffusionICfun1(x, P,kp_ref,D01)
%pdex1ic1 permet d'implémenter une condition initiale à l'equation partielle différentielle de l'éthylène
        %[M1(0); M2(0); T(0); Rpol(0)]
    u0 = [0;     0;     P.Tb; P.r_pol0];

    
end