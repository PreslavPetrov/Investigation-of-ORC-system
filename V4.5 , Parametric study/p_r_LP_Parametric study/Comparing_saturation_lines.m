clear
clc

out_water = Saturation_lines('Water');
out_R245fa = Saturation_lines('R245fa');

%plot(out_water(1,:),out_water(2,:),'--',out_R245fa(1,:),out_R245fa(2,:))
%title('Graph of s_s_a_t agaisntt t_s_a_t')                                    % Title
%xlabel('s_s_a_t [kJ/KgK]')                                                    % x-axis label
%ylabel('t_s_a_t [Degrees C]')                                                 % y-axis label
%legend('Water' ,'R245fa','Location','southeast')                              % Legend

T_cond_o = 20  ;

p_h = 20*(10^5);
t_sat_b_h = CoolProp.PropsSI('T' ,'P', p_h , 'Q' , 0, 'R245fa') - 273.15;

p_l = 2.5*(10^5);
t_sat_b_l = CoolProp.PropsSI('T' ,'P', p_l, 'Q' , 0, 'R245fa') - 273.15;

plot(out_R245fa(1,:),out_R245fa(2,:),'-',out_R245fa(1,:),t_sat_b_l,'*',out_R245fa(1,:),t_sat_b_h,'*',out_R245fa(1,:),T_cond_o,'o')
title('Graph of s_s_a_t agaisntt t_s_a_t')                                    % Title
xlabel('s_s_a_t [kJ/KgK]')                                                    % x-axis label
ylabel('t_s_a_t [Degrees C]')                                                 % y-axis label
legend('R245fa','Location','southeast')                                       % Legend



