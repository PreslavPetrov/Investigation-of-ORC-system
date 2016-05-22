clc
%POSTPROCESSING
s1 = size(out_ORC_opt);
s2 = size(xx);

%1.t_to_he_i  (4)
for j =1:1:s1(2)
for i = 1:1:s1(3)
ply(j,i) = out_ORC_opt(4,j,i);
end 
end
figure()
plot (xx,ply(1,:),'-r',xx,ply(2,:),'-g',xx,ply(3,:),'-b',xx,ply(4,:),'--');
grid on
title('t_t_ohe_i against various low cycle pressures of ORC System');       % Title
legend('50% Load', '75% Load','85% Load','100% Load');                      % Legend
xlabel('Low cycle pressure of ORC system [bar]');                           % x-axis label
ylabel('t_t_ohe_i [Degrees C] ');                                           % y-axis label

%2.t_sw_cond_o (8)
for j =1:1:s1(2)
for i = 1:1:s1(3)
ply2(j,i) = out_ORC_opt(8,j,i);
end 
end
figure()
plot (xx,ply2(1,:),'-r',xx,ply2(2,:),'-g',xx,ply2(3,:),'-b',xx,ply2(4,:),'--');
grid on
title('t_s_wcond_o against various low cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Low cycle pressure of ORC system [bar]');                            % x-axis label
ylabel('t_s_wcond_o [Degrees C] ');                                          % y-axis label

%3.t_g_o (15)
for j =1:1:s1(2)
for i = 1:1:s1(3)
ply3(j,i) = out_ORC_opt(15,j,i);
end 
end
figure()
plot (xx,ply3(1,:),'-r',xx,ply3(2,:),'-g',xx,ply3(3,:),'-b',xx,ply3(4,:),'--');
grid on
title('t_e_go against various low cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Low cycle pressure of ORC system [bar]');                            % x-axis label
ylabel('t_e_go [Degrees C] ');                                          % y-axis label

%4.t_to_he_o (16)
for j =1:1:s1(2)
for i = 1:1:s1(3)
ply4(j,i) = out_ORC_opt(16,j,i);
end 
end
figure()
plot (xx,ply4(1,:),'-r',xx,ply4(2,:),'-g',xx,ply4(3,:),'-b',xx,ply4(4,:),'--');
grid on
title('t_t_ohe_o against various low cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Low cycle pressure of ORC system [bar]');                            % x-axis label
ylabel('t_t_ohe_o [Degrees C] ');                                          % y-axis label

%5.mfr_to_b (10)
for j =1:1:s1(2)
for i = 1:1:s1(3)
ply5(j,i) = out_ORC_opt(10,j,i);
end 
end
figure()
plot (xx,ply5(1,:),'-r',xx,ply5(2,:),'-g',xx,ply5(3,:),'-b',xx,ply5(4,:),'--');
grid on
title('mfr_t_obp against various low cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Low cycle pressure of ORC system [bar]');                            % x-axis label
ylabel('mfr_t_obp [kg/s] ');                                              % y-axis label

%6.mfr_to_ev (11)
for j =1:1:s1(2)
for i = 1:1:s1(3)
ply6(j,i) = out_ORC_opt(11,j,i);
end 
end
figure()
plot (xx,ply6(1,:),'-r',xx,ply6(2,:),'-g',xx,ply6(3,:),'-b',xx,ply6(4,:),'--');
grid on
title('mfr_t_o_ev against various low cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Low cycle pressure of ORC system [bar]');                            % x-axis label
ylabel('mfr_t_oev [kg/s] ');                                              % y-axis label

%7.Q_he (17)
for j =1:1:s1(2)
for i = 1:1:s1(3)
ply7(j,i) = out_ORC_opt(17,j,i);
end 
end
figure()
plot (xx,ply7(1,:),'-r',xx,ply7(2,:),'-g',xx,ply7(3,:),'-b',xx,ply7(4,:),'--');
grid on
title('Q_h_e against various low cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Low cycle pressure of ORC system [bar]');                            % x-axis label
ylabel('Q_h_e, Exhaust gas heat exchanger thermal Power [kW] ');                                              % y-axis label

%8.Q_ev (5)
for j =1:1:s1(2)
for i = 1:1:s1(3)
ply8(j,i) = out_ORC_opt(5,j,i);
end 
end
figure()
plot (xx,ply8(1,:),'-r',xx,ply8(2,:),'-g',xx,ply8(3,:),'-b',xx,ply8(4,:),'--');
grid on
title('Q_e_v  against various low cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Low cycle pressure of ORC system [bar]');                            % x-axis label
ylabel('Q_e_v, Evaporator thermal Power [kW]  ');                                              % y-axis label

%9.P_el (12)
for j =1:1:s1(2)
for i = 1:1:s1(3)
ply9(j,i) = out_ORC_opt(12,j,i);
end 
end
figure()
plot (xx,ply9(1,:),'-r',xx,ply9(2,:),'-g',xx,ply9(3,:),'-b',xx,ply9(4,:),'--');
grid on
title('P_e_l  against various low cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Low cycle pressure of ORC system [bar]');                            % x-axis label
ylabel('P_e_l, Generator produced electrical power [kW] ');                                              % y-axis label

%10.P_net (13)
for j =1:1:s1(2)
for i = 1:1:s1(3)
ply10(j,i) = out_ORC_opt(13,j,i);
end 
end
figure()
plot (xx,ply10(1,:),'-r',xx,ply10(2,:),'-g',xx,ply10(3,:),'-b',xx,ply10(4,:),'--');
grid on
title('P_n_e_t  against various low cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Low cycle pressure of ORC system [bar]');                            % x-axis label
ylabel('P_n_e_t , Net Electrical Power [kW] ');                                              % y-axis label

%11.eff_turb_is (6) 
for j =1:1:s1(2)
for i = 1:1:s1(3)
ply10(j,i) = out_ORC_opt(6,j,i);
end 
end
figure()
plot (xx,ply10(1,:),'-r',xx,ply10(2,:),'-g',xx,ply10(3,:),'-b',xx,ply10(4,:),'--');
grid on
title('eff_t_u_r_bis  against various low cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Low cycle pressure of ORC system [bar]');                            % x-axis label
ylabel('eff_t_u_r_bis, Isentropic turbine efficiency [-] ');                                              % y-axis label

%12.eff_turb_m (14)
for j =1:1:s1(2)
for i = 1:1:s1(3)
ply10(j,i) = out_ORC_opt(14,j,i);
end 
end
figure()
plot (xx,ply10(1,:),'-r',xx,ply10(2,:),'-g',xx,ply10(3,:),'-b',xx,ply10(4,:),'--');
grid on
title('eff_t_u_r_bm  against various low cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Low cycle pressure of ORC system [bar]');                            % x-axis label
ylabel('eff_t_u_r_bm , Mechanical Turbine-Generator efficiency [-] ');                                              % y-axis label