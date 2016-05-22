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
vfr_to = 0.01605498965932;     %m^3 /s
vfr_sw = 0.052538324106033;    %m^3 /s
vfr_r = 0.007698867314537;     %m^3/s    
Refrigerant = 'R245fa' 

    % Optimizing    
t_to_he_i = [120, 120, 120, 120];                                                             % Temperature of Thermal oil at heat exchanger inlet [Degrees C]
UA_he_d = 30.4975;                                                                % UA of exhasut gas Heat exchanger , UA = NTU*Cmin [kW/m^2]
UA_ev_d = 31.9304;

%Run supportign databases
thermal_properties_data_to
eff_pump_data

%1. Exhaust gas heat exchanger
t_g_o_est = 180;                                                             % Initial guess of Temperature of Exhasut gas at outlet, only used for estimation of CP for the exhaust gas [Degrees C]
t_to_he_o_est = 300;                                                         % Initial guess of Temperature of Thermal oil at outlet, only used for estimation of CP for the Thermal oil [Degrees C]
t_g_min_req = 0.1239e6*s^3 - 1.3644e4*s^2 + 6.5509e2*s + 124.97;             % Minimum Allowable temperature 


for j =1:1:4
t_g_o_est = 180;
 for i =1:1:5
    k(j,:) = fun_compos_g(c, h, s, lamda(j));                                % composition of the exhaust gas (kg i/kg of exhaust gas): k = [k_CO2, k_H2O, k_SO2, k_O2, k_N2]
% Estimation of cpm_g  
    cpm_g(j) = fun_cpm_g(t_g_o_est,t_g_i(j),k(j,:));                         % Estiamtion of Mean specific heat at constant pressure of exhaust gas [kJ/kgK]
    C_g(j)= mfr_g(j)*cpm_g(j);                                               % Estimation of Exhast gas Heat capacity [KJ/sK = kW/K]
        
% Estimation of Cp_to_he     
    Cp_to_ev = pchip(TO_Properties(:,1),TO_Properties(:,2), (t_to_he_o_est+t_to_he_i(j))/2); % Estiamtion of Mean specific heat at constant pressure of Thermal oil in Exhasut gas heat exchanger [kJ/kgK]
    Dens_to_ev(j) = pchip(TO_Properties(:,1),TO_Properties(:,3), (t_to_he_o_est+t_to_he_i(j))/2);
    mfr_to(j) = vfr_to*Dens_to_ev(j);
    C_to_he(j) = mfr_to(j)*Cp_to_ev;                                            % Estimation of Thermal oil Heat capacity at exhasut gas heat exchanger [KJ/sK = kW/K]

% Correction for UA for off-design operation 
    UA_he_p(j) = UA_he_d*(mfr_g(j)/mfr_g(3))^0.65;

% From using MLTD method and Energy Balance:      
Co = C_g/C_to_he;


F(1) = log((t_g_i(j)-t_to_he_o(j))/(t_g_o_est-t_to_he_i(j))) - UA_he_p*(t_g_i(j)- t_to_he_o(j))/(C_g*(t_g_i(j)-t_g_o_est)) ;  
F(2) = Co*(t_g_i(j) - t_g_o(j)) + t_to_he_i(j) - t_to_he_o(j)
 
 end
 
end   
t_g_o
t_to_he_o
%out_ORC


