clc
clear
%INPUT
   % Design
c = 0.86;                                                                    % Carbon mass fraction [kg C/kg fuel] 
h = 0.11;                                                                    % Hudrogen mass fraction [kg H/kg fuel]
s = 0.03;                                                                    % Sulfur mass fraction [kg S/kg fuel]
LHV = 42700;                                                                 % Fuel lowe heating value [kJ/kg]
lamda = [2.172223221 2.215813297 2.311768258 2.215565168];                   % Exhaust gas air to fuel equivalence ratio [-]
mfr_g = [7.566	11.252	13.289	15.423*0.8];                                 % Exhaust gas mass flow rate [kg/s]
t_g= [415; 350; 330; 390];                                                   % Temperature of exhaust gas exiting the engine exhaust gas turbochacher  [Degrees C]
Dt_g = 4 ;                                                                   % Temperature drop in the exhaust pipe (connecting the engine t/c turbine and the boiler (usually <= 5 deg C)
t_g_i = t_g - Dt_g;                                                          % Temperature of exhasut gas entering the exhaust gas heat exchanger [Degrees C]
vfr_to_pump = 0.015473488755932;     %m^3/s
vfr_sw_pump = 0.052362903188870;    %m^3/s
vfr_r_pump = 0.007698867314537;     %m^3/s     
Refrigerant = 'R245fa' 
L_tg_rated = 400;   % Rated Genrator power [kW]                                                         % Rate Gen load [kW]

    % Parametrs that could be optimized  
t_to_he_i = [120, 120, 120, 120];                                           % Temperature of Thermal oil at heat exchanger inlet [Degrees C]
UA_he_d = 30.4975;                                                                % UA of exhasut gas Heat exchanger , UA = NTU*Cmin [kW/m^2]
UA_ev_d = 31.9283;
p_r_LP_bar = 2.5;                                                           % Lower Pressure of ORC system [bar]
p_r_UP_bar = 20;                                                            % Upper Pressure of ORC system [bar]
  
% Optimized Parameter
p_r_UP_bar_min = 20;                                                           % Lower Pressure of ORC system [bar]
l = 0;
p_r_UP_bar_max  = CoolProp.PropsSI('PCRIT' ,'P', p_r_LP_bar*(10^5) , 'Q' ,1, 'R245FA')/(10^5);

for p_r_UP_bar = (p_r_LP_bar+1.5):1:p_r_UP_bar_max
l = l+1;
xx(l) = p_r_UP_bar;
t_to_he_i = [120, 120, 120, 120];  
%Run supportign databases
thermal_properties_data_to
eff_pump_data

%1. Exhaust gas heat exchanger
t_g_min_req = 0.1239e6*s^3 - 1.3644e4*s^2 + 6.5509e2*s + 124.97;             % Minimum Allowable temperature 
p_r_LP = p_r_LP_bar*100000;                                                 % Lower Pressure of ORC system [Pa], needed for CoolProp
p_r_UP = p_r_UP_bar*100000;                                                 % Upper Pressure of ORC system [Pa], needed for CoolProp
out_ORC = zeros(14,4);
for j =1:1:4
    t_to_he_o_est = 150;                                                         % Initial guess of Temperature of Thermal oil at outlet, only used for estimation of CP for the Thermal oil [Degrees C]
    t_g_o_est = 180;                                                         % Initial guess of Temperature of Exhasut gas at outlet, only used for estimation of CP for the exhaust gas [Degrees C]
 for i =1:1:20
    k(j,:) = fun_compos_g(c, h, s, lamda(j));                                % composition of the exhaust gas (kg i/kg of exhaust gas): k = [k_CO2, k_H2O, k_SO2, k_O2, k_N2]
% Estimation of cpm_g  
    cpm_g(j) = fun_cpm_g(t_g_o_est,t_g_i(j),k(j,:));                         % Estiamtion of Mean specific heat at constant pressure of exhaust gas [kJ/kgK]
    C_g(j)= mfr_g(j)*cpm_g(j);                                               % Estimation of Exhast gas Heat capacity [KJ/sK = kW/K]
        
% Estimation of massflworate of themral oil.
    Dens_to_pump = pchip(TO_Properties(:,1),TO_Properties(:,3), t_to_he_i(j));
    mfr_to(j) = vfr_to_pump*Dens_to_pump;

% Estimation of Cp_to_he     
    Cp_to_he = pchip(TO_Properties(:,1),TO_Properties(:,2), (t_to_he_o_est+t_to_he_i(j))/2); % Estiamtion of Mean specific heat at constant pressure of Thermal oil in Exhasut gas heat exchanger [kJ/kgK]  
    C_to_he(j) = mfr_to(j)*Cp_to_he;                                            % Estimation of Thermal oil Heat capacity at exhasut gas heat exchanger [KJ/sK = kW/K]

    Co(j) = C_g(j)/C_to_he(j);
% Correction for UA for off-design operation 
    UA_he_p(j) = UA_he_d*(mfr_g(j)/mfr_g(3))^0.65;
csvwrite('UA_he_p.txt',UA_he_p(j))
csvwrite('t_to_he_i.txt',t_to_he_i(j))
csvwrite('t_g_i.txt',t_g_i(j))
csvwrite('C_g.txt',C_g(j))
csvwrite('Co.txt',Co(j))

% From using MLTD method and Energy Balance:
fun = @fun_HE_MLTD;
x0 = [0,0];
x = fsolve(fun,x0);
t_to_he_o(j) = real(x(1));
t_g_o(j) = real(x(2));
            
    t_to_he_o_est = t_to_he_o(j) ;                                             % Assigning the new t_to_he_o value to t_to_he_o_guess
    t_g_o_est = t_g_o(j) ;                                                      % Assigning the new t_g_o value to  t_g_o_est

 
Q_he(j) = C_g(j)*(t_g_i(j)-t_g_o(j));
Dp_to_pump_bar(j) = 1.728-1.2+ 0.5*Dens_to_pump*9.81/(10^5);                 % Pressure difference that the thermal oil pump needs to provide   //Lowsest pressure in the system is chosen as 1.2 bar, hence if 20% losses of pressure are assuemd in the heaty exchanegr and evaporator this leads to 1.728 bar for the upper pressure i nthe system. Since manifacturers of pumps recommend an additional 0.5 m (~0.044145 bar) NPSHR margins                                
eff_to_pump = pchip(eff_max_pump_data(:,1),eff_max_pump_data(:,2), vfr_to_pump);  % Estimation of efficiency of thermal oil pump
P_to_pump(j) = vfr_to_pump*Dp_to_pump_bar(j)*(10^2)/(eff_to_pump);                    % Power of Thermal oil pump (kW) // W = vfr*Dp

out_ORC(:,j) = fun_ORC(t_to_he_o(j),p_r_LP_bar,p_r_UP_bar,vfr_r_pump,UA_ev_d,mfr_to(j),t_to_he_i(j),Refrigerant,vfr_sw_pump,P_to_pump(j),L_tg_rated);
t_to_he_i(j) = out_ORC (4,j);
mfr_to_b(j) = out_ORC (10,j);
mfr_to_ev(j) = out_ORC (11,j);
 end
 %CHECKS
 % Sulphuric acid formation check
 if t_g_o(j) < t_g_min_req 
     disp(['The temperature of the exhasut gas, ' num2str(t_g_o(j)) ', at the outlet is lower than ', num2str(t_g_min_req ) , ' deg C and should be increased in order to avoid sulphuric acid formation'])
 end
 
%mfr_b/mfr_to check
 if mfr_to_b/mfr_to >= 1 
     disp([' Bled massflow rate should be in the range of (0,1)'])
 end
 if mfr_to_b/mfr_to <= 0
     disp([' Bled massflow rate should be in the range of (0,1)'])
 end
end
out_ORC(15,:) = t_g_o;
out_ORC(16,:) = t_to_he_o;
out_ORC(17,:) = Q_he;

VH = out_ORC(6,:);
VR = out_ORC(7,:);
%eff_m = out_ORC(14,:);
out_ORC_opt(:,:,l) = out_ORC(:,:); 

end

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
title('t_t_ohe_i against various upper cycle pressures of ORC System');       % Title
legend('50% Load', '75% Load','85% Load','100% Load');                      % Legend
xlabel('Upper cycle pressure of ORC system [bar]');                           % x-axis label
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
title('t_s_wcond_o against various upper cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Upper cycle pressure of ORC system [bar]');                            % x-axis label
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
title('t_e_go against various upper cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Upper cycle pressure of ORC system [bar]');                            % x-axis label
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
title('t_t_ohe_o against various upper cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Upper cycle pressure of ORC system [bar]');                            % x-axis label
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
title('mfr_t_obp against various upper cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Upper cycle pressure of ORC system [bar]');                            % x-axis label
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
title('mfr_t_o_ev against various upper cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Upper cycle pressure of ORC system [bar]');                            % x-axis label
ylabel('mfr_t_o_ev [kg/s] ');                                              % y-axis label

%7.Q_he (17)
for j =1:1:s1(2)
for i = 1:1:s1(3)
ply7(j,i) = out_ORC_opt(17,j,i);
end 
end
figure()
plot (xx,ply7(1,:),'-r',xx,ply7(2,:),'-g',xx,ply7(3,:),'-b',xx,ply7(4,:),'--');
grid on
title('Q_h_e against various upper cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Upper cycle pressure of ORC system [bar]');                            % x-axis label
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
title('Q_e_v  against various upper cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Upper cycle pressure of ORC system [bar]');                            % x-axis label
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
title('P_e_l  against various upper cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Upper cycle pressure of ORC system [bar]');                            % x-axis label
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
title('P_n_e_t  against various upper cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Upper cycle pressure of ORC system [bar]');                            % x-axis label
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
title('eff_t_u_r_bis  against various upper cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Upper cycle pressure of ORC system [bar]');                            % x-axis label
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
title('eff_t_u_r_bm  against various upper cycle pressures of ORC System');      % Title
legend('50% Load', '75% Load','85% Load','100% Load');                       % Legend
xlabel('Upper cycle pressure of ORC system [bar]');                            % x-axis label
ylabel('eff_t_u_r_bm , Mechanical Turbine-Generator efficiency [-] ');                                              % y-axis label

csvwrite('out_ORC_opt2.txt',out_ORC_opt)
