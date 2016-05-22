function [eff_mech] = fun_mech_eff_turb (L_tg_rated,El_load_est)
Gen_mech_eff_data
%input
%L_tg_rated = 280
%El_load_est = 80

% eff_mech_rated
if  L_tg_rated < 290 
    x = gen_eff_m (:,1);
    y = gen_eff_m (:,2);
    n = 10;
    Poly_mech_eff_rated=polyfit(x,y,n);

    eff_mech_rated = 0 ; 
        for i =1:1:n+1
         eff_mech_rated = eff_mech_rated + Poly_mech_eff_rated(i)*(L_tg_rated)^(n-i+1);
        end
else
    eff_mech_rated = pchip(gen_eff_m(:,1),gen_eff_m(:,2), L_tg_rated);
end
% Load_cor_f
lcf = pchip(Load_cor_f(:,1),Load_cor_f(:,2), El_load_est);

% eff_mech
eff_mech = eff_mech_rated * lcf;

end