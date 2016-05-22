function F  = fun_HE_MLTD(x)
%x(1) t_to_he_o
%x(2) t_g_he_o

UA_he_p = csvread('UA_he_p.txt');
t_to_he_i = csvread('t_to_he_i.txt');
t_g_i = csvread('t_g_i.txt');
C_g = csvread('C_g.txt');
Co = csvread('Co.txt');


F(1)  = log((t_g_i-x(1))/(x(2)-t_to_he_i)) - UA_he_p*((t_g_i - x(1))-(x(2)-t_to_he_i))/(C_g*(t_g_i-x(2)));  
F(2) = Co*(t_g_i - x(2)) + t_to_he_i - x(1);

end
