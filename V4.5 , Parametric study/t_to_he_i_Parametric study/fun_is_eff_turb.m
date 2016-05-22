function [eff_is] = fun_is_eff_turb (Ns)
Gen_is_eff_data
%input
%Ns

%eff_is

if  Ns < 0.25 
    x = gen_eff_is (:,1);
    y = gen_eff_is (:,2);
    n = 1;
    Poly_is_eff_rated=polyfit(x,y,n);

    eff_is = 0 ; 
        for i =1:1:n+1
         eff_is = eff_is + Poly_is_eff_rated(i)*(Ns)^(n-i+1);
        end
else
    eff_is = pchip(gen_eff_is(:,1),gen_eff_is(:,2), Ns);
end


end