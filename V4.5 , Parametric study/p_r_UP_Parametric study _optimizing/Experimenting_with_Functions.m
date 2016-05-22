clear
clc

p_r_LP_bar = 1;
T  = CoolProp.PropsSI('T' ,'P', p_r_LP_bar*(10^5), 'Q' ,0, 'R245FA') - 273.15