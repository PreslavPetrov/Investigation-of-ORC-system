function [out] = Saturation_lines(Refrigerant)
%Refrigerant = 'Water'
p_crit =  CoolProp.PropsSI('PCRIT' ,'P', 100000, 'Q' , 1, Refrigerant);
p_LP_bar =1 ;
p_HP_bar =p_crit/100000;
p_LP = p_LP_bar*100000;
p_HP = p_HP_bar*100000;

pd= 100000;
i = 0;
for pi = p_LP:pd:p_HP
i =i+1;
t_sat(i)= CoolProp.PropsSI('T' ,'P', pi, 'Q' , 0, Refrigerant)-273.15;
s_sat(i)= CoolProp.PropsSI('S' ,'P', pi, 'Q' , 0, Refrigerant)/(10^3);
end

for pi = p_HP:-pd:p_LP
i =i+1;
t_sat(i)= CoolProp.PropsSI('T' ,'P', pi, 'Q' , 1, Refrigerant) -273.15;
s_sat(i)= CoolProp.PropsSI('S' ,'P', pi, 'Q' , 1, Refrigerant)/(10^3);
end

%plot(s_sat,t_sat,'--')
%title('Graph of s_s_a_t agaisntt t_s_a_t')                                    % Title
%xlabel('s_s_a_t [kJ/KgK]')                                                    % x-axis label
%ylabel('t_s_a_t [Degrees C]')                                                % y-axis label
%legend('Saturation line','Location','southeast')                              % Legend

out(1,:) = s_sat;
out(2,:) = t_sat;
end
