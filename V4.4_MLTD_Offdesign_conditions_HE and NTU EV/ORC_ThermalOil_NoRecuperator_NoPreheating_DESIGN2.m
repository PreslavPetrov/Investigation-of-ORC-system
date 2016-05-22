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
L_tg_rated = 400;                                                           % Rate Gen load [kW]

    % Optimzing    
t_to_he_i = [120, 120, 120, 120];                                           % Temperature of Thermal oil at heat exchanger inlet [Degrees C]
UA_he_d = 30.4975;                                                                % UA of exhasut gas Heat exchanger , UA = NTU*Cmin [kW/m^2]
UA_ev_d = 31.9283;
p_r_LP_bar = 2.5;                                                           % Lower Pressure of ORC system [bar]
p_r_UP_bar = 20;                                                            % Upper Pressure of ORC system [bar]
p_r_LP = p_r_LP_bar*100000;                                                 % Lower Pressure of ORC system [Pa], needed for CoolProp
p_r_UP = p_r_UP_bar*100000;                                                 % Upper Pressure of ORC system [Pa], needed for CoolProp

%Run supportign databases
thermal_properties_data_to
eff_pump_data

%1. Exhaust gas heat exchanger
t_g_min_req = 0.1239e6*s^3 - 1.3644e4*s^2 + 6.5509e2*s + 124.97;             % Minimum Allowable temperature 
 
for j =1:1:4
    t_to_he_o_est = 140;                                                         % Initial guess of Temperature of Thermal oil at outlet, only used for estimation of CP for the Thermal oil [Degrees C]
    t_g_o_est = 200;                                                         % Initial guess of Temperature of Exhasut gas at outlet, only used for estimation of CP for the exhaust gas [Degrees C]
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
t_to_he_i(j) = out_ORC(4,j);
mfr_to_b = out_ORC(8,:);
mfr_to_ev = out_ORC(9,j);
 end
 %CHECKS
 % Sulphuric acid formation check
 if t_g_o(j) < t_g_min_req 
     disp(['The temperature of the exhasut gas, ' num2str(t_g_o(j)) ', at the outlet is lower than ', num2str(t_g_min_req ) , ' deg C and should be increased in order to avoid sulphuric acid formation'])
 end
 
%mfr_b/mfr_to check
 if mfr_to_b/mfr_to >= 1 
     disp([' Bled massflow rate should be in the range of (0,1)'])
     mfr_to_b/mfr_to
 end
 if mfr_to_b/mfr_to <= 0
     disp([' Bled massflow rate should be in the range of (0,1)'])
     mfr_to_b/mfr_to
 end
end
out_ORC(14,:) = t_g_o;
out_ORC(15,:) = t_to_he_o;
out_ORC

