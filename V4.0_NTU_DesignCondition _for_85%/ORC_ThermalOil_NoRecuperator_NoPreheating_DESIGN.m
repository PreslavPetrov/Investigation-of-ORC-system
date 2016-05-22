clc
clear
%INPUT
c = 0.86;                                                                    % Carbon mass fraction [kg C/kg fuel] 
h = 0.11;                                                                    % Hudrogen mass fraction [kg H/kg fuel]
s = 0.03;                                                                    % Sulfur mass fraction [kg S/kg fuel]
LHV = 42700;                                                                 % Fuel lowe heating value [kJ/kg]
lamda = [2.172223221 2.215813297 2.311768258 2.215565168];                   % Exhaust gas air to fuel equivalence ratio [-]
mfr_g = [7.566	11.252	13.289	15.423*0.8];                                 % Exhaust gas mass flow rate [kg/s]
t_g= [415; 350; 330; 390];                                                   % Temperature of exhaust gas exiting the engine exhaust gas turbochacher  [Degrees C]
Dt_g = 4 ;                                                                   % Temperature drop in the exhaust pipe (connecting the engine t/c turbine and the boiler (usually <= 5 deg C)
t_g_i = t_g - Dt_g;                                                          % Temperature of exhasut gas entering the exhaust gas heat exchanger [Degrees C]

t_to_he_i = 120;                                                             % Temperature of Thermal oil at heat exchanger inlet [Degrees C]                                                               % UA of exhasut gas Heat exchanger , UA = NTU*Cmin [kW/m^2]
eff_he = 0.80;
%Dt_superheat_ev = 10 ;                                                       % Degree of Superheat of Refrigerant at Evaporator outlet [-]
Refrigerant = 'R245fa' 
mfr_to = 15;

%Run supportign databases
eff_cor_factors_data
thermal_properties_data_to
eff_pump_data

%1. Exhaust gas heat exchanger
t_g_o_est= 120;                                                              % Initial guess of Temperature of Exhasut gas at outlet, only used for estimation of CP for the exhaust gas [Degrees C]
t_to_he_o_est = 300;                                                       % Initial guess of Temperature of Thermal oil at outlet, only used for estimation of CP for the Thermal oil [Degrees C]
UA = 0;

% ITERATING FOR ESTIAMTION OF CpS
    for i =1:1:3        
    k = fun_compos_g(c, h, s, lamda(3));                    % composition of the exhaust gas (kg i/kg of exhaust gas): k = [k_CO2, k_H2O, k_SO2, k_O2, k_N2]
% Estimation of cpm_g  
    cpm_g = fun_cpm_g(t_g_o_est,t_g_i(3),k);                         % Estiamtion of Mean specific heat at constant pressure of exhaust gas [kJ/kgK]
    C_g= mfr_g(3)*cpm_g;                                              % Estimation of Exhast gas Heat capacity [KJ/sK = kW/K]
        
% Estimation of Cp_to_he     
    Cp_to_ev = pchip(TO_Properties(:,1),TO_Properties(:,2), (t_to_he_o_est+t_to_he_i)/2); % Estiamtion of Mean specific heat at constant pressure of Thermal oil in Exhasut gas heat exchanger [kJ/kgK]
    C_to_he = mfr_to*Cp_to_ev;                                               % Estimation of Thermal oil Heat capacity at exhasut gas heat exchanger [KJ/sK = kW/K]

% From using Effectiveness-NTU method and Energy Balance:
         if(C_g>C_to_he)   % Cmin = C_to_he
             t_g_o = t_g_i(3) - (t_g_i(3)-t_to_he_i)*eff_he*C_to_he/C_g;
             t_to_he_o = t_to_he_i + C_g*(t_g_i(3)-t_g_o)/C_to_he;
             
             Cr = C_to_he/C_g;                                         % Heat capacity ratio [-], as Cr decreases , effectiveness increases
             NTU_he = (1/(1-Cr))*log((1-eff_he*(Cr))/(1-eff_he));
             UA_he =NTU_he*C_to_he;
         else   % Cmin = C_g       
             t_g_o = t_g_i(3) - (t_g_i(3)-t_to_he_i)*eff_he;
             t_to_he_o = t_to_he_i +  C_g*(t_g_i(3)-t_g_o)/C_to_he; 
             
             Cr = C_g/C_to_he ;                                            
             NTU_he = (1/(1-Cr))*log((1-eff_he*(Cr))/(1-eff_he));
             UA_he =NTU_he*C_g;
         end
    t_to_he_o_est = t_to_he_o ;                                            % Assigning the new t_to_he_o value to t_to_he_o_guess
    t_g_o_est = t_g_o;                                                       % Assigning the new t_g_o value to  t_g_o_est
    Q_he = C_g*(t_g_i(3)-t_g_o);

   end
                                                                             % Exhasut gas heat exchanger Power (kW)
den_to_he_i = pchip(TO_Properties(:,1),TO_Properties(:,3), t_to_he_i);       % Estiamtion of Density at constant pressure of Thermal oil passing throguh themral oil pump [kg/m^3]
Dp_to_pump_bar = 1.728-1.2+ 0.5*den_to_he_i*9.81/(10^5);                     % Pressure difference that the thermal oil pump needs to provide   //Lowsest pressure in the system is chosen as 1.2 bar, hence if 20% losses of pressure are assuemd in the heaty exchanegr and evaporator this leads to 1.728 bar for the upper pressure i nthe system. Since manifacturers of pumps recommend an additional 0.5 m (~0.044145 bar) NPSHR margins
0.5*den_to_he_i*9.81/(10^5)
vfr_to_pump = mfr_to/den_to_he_i ;                                        % Volume flow rate of Thermal oil through Themral oil pump[m^3 /s]
eff_to_pump = pchip(eff_max_pump_data(:,1),eff_max_pump_data(:,2), vfr_to_pump);  % Estimation of efficiency of thermal oil pump
P_to_pump = vfr_to_pump*Dp_to_pump_bar*(10^2)/(eff_to_pump) ;       % Power of Thermal oil pump (kW) // W = vfr*Dp

%2. Interm Loop
t_to_ev_i = t_to_he_o;                                                       % Inlet Temperature of Thermal oil to Evaporator [Degrees C]

%4. ORC
%PUMP
p_r_LP_bar = 2.5;                                                                % Lower Pressure of ORC system [bar]
p_r_UP_bar = 20;                                                             % Upper Pressure of ORC system [bar]
p_r_LP = p_r_LP_bar*100000;                                                      % Lower Pressure of ORC system [Pa], needed for CoolProp
p_r_UP = p_r_UP_bar*100000;                                                  % Upper Pressure of ORC system [Pa], needed for CoolProp


tK_r_cond_sat = CoolProp.PropsSI('T' ,'P', p_r_LP, 'Q' , 0, 'R245FA');       % Saturation temeprature of Refrigerant at Lower pressure, p_r_LP [K]
t_r_cond_sat = tK_r_cond_sat - 273.15                                       % Saturation temeprature of Refrigerant at Lower pressure, p_r_LP [Degrees C]
t_r_cond_o = t_r_cond_sat;                                                    % Outlet Temperature of Refrigerant at Condenser [Degrees C]
tK_r_cond_o = t_r_cond_o + 273.15;                                           % Outlet Temperature of Refrigerant at Condenser [K], needed for CoolProp
s_r_cond_o = CoolProp.PropsSI('S' ,'P', p_r_LP, 'Q' , 0, 'R245FA');          % Outlet Entropy of Refrigerant at Condneser [J/kgK]
h_r_cond_o = CoolProp.PropsSI('H' ,'P', p_r_LP, 'Q' , 0, 'R245FA');          % Outlet Enthalpy of Refrigerant at Condneser [J/kg]
s_r_pump_o_is = s_r_cond_o;                                                  % Isentropic Outlet Entropy of Refrigerant at Refrigerant Pump [J/kgK]
h_r_pump_o_is = CoolProp.PropsSI('H' ,'P', p_r_UP, 'S' , s_r_pump_o_is, 'R245FA');  % Isentropic Outlet Enthalpy of Refrigerant at Refrigerant Pump [J/kg]
eff_is_pump_r = 0.80;                                                        % Isentropic efficiency of the Refrigerant pump [-]
h_r_pump_o = h_r_cond_o + (h_r_pump_o_is-h_r_cond_o)/eff_is_pump_r;          % Outlet Enthalpy of Refrigerant at Refrigerant Pump [J/kg]
s_r_pump_o = CoolProp.PropsSI('S' ,'P', p_r_UP, 'H' , h_r_pump_o, 'R245FA');

tK_r_pump_o = CoolProp.PropsSI('T' ,'P', p_r_UP, 'H' , h_r_pump_o, 'R245FA'); % Outlet Temperature of Refrigerant at Refrigerant Pump [K]
t_r_pump_o = tK_r_pump_o -273.15; 
dens_r_pump_o = CoolProp.PropsSI('D' ,'P', p_r_UP, 'H' , h_r_pump_o, 'R245FA');% Outlet Density of Refrigerant at Refrigerant Pump [K]
dens_r_pump_i = CoolProp.PropsSI('D' ,'P', p_r_LP, 'H' , h_r_cond_o, 'R245FA');% Inlet Density of Refrigerant at Refrigerant Pump [K]
dens_r_pump = (dens_r_pump_o + dens_r_pump_i)/2;                             % Average Density of Refrigerant at Refrigerant Pump [K], note that Dt in pump is about 1 degree

mfr_r = 10; 
vfr_r_pump = mfr_r/dens_r_pump ;                                             % Requred Volume flowrate of Refrigerant [m^3/s]
eff_r_pump =  pchip(eff_max_pump_data(:,1),eff_max_pump_data(:,2), vfr_r_pump);  % Estimation of efficiency of refrigerant pump [-]
P_pump_r = vfr_r_pump*(p_r_UP-p_r_LP)/(eff_r_pump*(10^3));   


%EVAPORATOR
eff_ev = 0.80;
t_to_ev_o_est = 80;
t_r_ev_i = t_r_pump_o;
h_r_ev_i = h_r_pump_o/(10^3);   %in kJ/kg
for i =1:1:5 
Cp_to_ev = pchip(TO_Properties(:,1),TO_Properties(:,2), (t_to_ev_o_est+t_to_ev_i)/2);
Cp_to_ev_min = pchip(TO_Properties(:,1),TO_Properties(:,2), (t_r_ev_i+t_to_ev_i)/2);
t_to_ev_o = t_to_ev_i - (Cp_to_ev_min/Cp_to_ev)*(t_to_ev_i-t_r_ev_i)*eff_ev;
t_to_ev_o_est = t_to_ev_o;

Cp_to_he_i = pchip(TO_Properties(:,1),TO_Properties(:,2), t_to_he_i);
Cp_to_he_o = pchip(TO_Properties(:,1),TO_Properties(:,2), t_to_he_o);
Cp_to_ev_o = pchip(TO_Properties(:,1),TO_Properties(:,2), t_to_ev_o);
mfr_to_ev = mfr_to*(Cp_to_he_i*t_to_he_i - Cp_to_he_o*t_to_he_o)/(Cp_to_ev_o*t_to_ev_o - Cp_to_he_o*t_to_he_o);
mfr_to ;
mfr_to_b = mfr_to - mfr_to_ev;
mfr_to_ev;
%mfr_to*Cp_to_he_i*t_to_he_i - mfr_to_b*Cp_to_he_o*t_to_he_o -mfr_to_ev*Cp_to_ev_o*t_to_ev_o 

h_r_ev_o = h_r_ev_i + mfr_to_ev*Cp_to_ev*(t_to_ev_i - t_to_ev_o)/mfr_r;  %kJ/kg

t_r_ev_o = CoolProp.PropsSI('T' ,'P', p_r_UP, 'H' , h_r_ev_o*(10^3), 'R245FA') -273.15;
h_r_ev_sat = CoolProp.PropsSI('H' ,'P', p_r_UP, 'Q' , 1, 'R245FA')/(10^3);
t_r_ev_sat = CoolProp.PropsSI('T' ,'P', p_r_UP, 'Q' , 1, 'R245FA') -273.15;
Dt_ev = t_r_ev_o - t_r_ev_sat;
end

if Dt_ev > 15
    'Degree of superheat should be between 0 and 15 degrees and it '
    Dt_ev 
    ' .Consider increasing the refrigerant mass  flowrate!'        
end

if h_r_ev_sat > h_r_ev_o
    'The Organic fluid ahs not evaporated compeltely in the evaporator! Consider reducing the mass flowrade of refrigerant'    
end 
Q_ev = mfr_to_ev*Cp_to_ev*(t_to_ev_i - t_to_ev_o);
NTU_ev = log(1/(1-eff_ev));
UA_ev = NTU_ev*Cp_to_ev*mfr_to_ev ;


%TURBINE
h_r_turb_i = h_r_ev_o;                                                       % Inlet Enthalpy of Refrigerant at Turbine inlet [J/kg]
t_r_turb_i = t_r_ev_o;
tK_r_turb_i = t_r_turb_i +273.15;
s_r_turb_i = CoolProp.PropsSI('S' ,'P', p_r_UP, 'T' , tK_r_turb_i, Refrigerant);  % Inlet Entropy of Refrigerant at Turbine [J/kgK]                                                  % 
s_r_turb_o_is = s_r_turb_i;                                                  % Isentropic Outlet Entropy of Refrigerant at Turbine [J/kgK]    
h_r_turb_o_is = CoolProp.PropsSI('H' ,'P', p_r_LP, 'S' , s_r_turb_o_is, Refrigerant); % Isentropic Outlet Enthalpy of Refrigerant at Turbine [J/kg]    

%ISENTROPIC EFFICIENCY OF TURBINE
N_turb = 20000 ; % RPM of turbine in rev per min
Ns_turb = N_turb/60 ; % RPS of turbine in rev per sec
Dens_r_turb_o_is = CoolProp.PropsSI('D' ,'P', p_r_LP, 'H' , h_r_turb_o_is, Refrigerant);
vfr_r_turb_o = mfr_r/(Dens_r_turb_o_is);
delh_is = (h_r_turb_i*(10^3)) - h_r_turb_o_is;
Ns = Ns_turb*((vfr_r_turb_o)^0.5)/(delh_is^(0.75));
eff_turb_is = fun_is_eff_turb(Ns) ;                                                             % Isentropic Efficiency of Turbine

h_r_turb_o = h_r_turb_i - eff_turb_is*(h_r_turb_i - h_r_turb_o_is/(10^3));          % Outlet Entropy of Refrigerant at Turbine [J/kgK]   
t_r_turb_o = CoolProp.PropsSI('T' ,'P', p_r_LP, 'H' , h_r_turb_o*(10^3), Refrigerant) -273.15;
s_r_turb_o = CoolProp.PropsSI('S' ,'P', p_r_LP, 'H' , h_r_turb_o*(10^3), Refrigerant);
              
L_tg_rated = 400;
El_load_est = 80;
for i =1:1:20
eff_m  = fun_mech_eff_turb (L_tg_rated,El_load_est);
P_el = mfr_r*(h_r_turb_i - h_r_turb_o)*eff_m ;                            % Isentropic Power of Turbine (REfrigerant power) [kJ/s = kW]
El_load_est =100*P_el/L_tg_rated;
end
P_el = mfr_r*(h_r_turb_i - h_r_turb_o)*eff_m ;                           % Power of Turbine (REfrigerant power) [kJ/s = kW]

%CONDENSER
Q_cond =  mfr_r*(h_r_turb_o - h_r_cond_o/(10^3)) ;                           % Condenser Thermodynamic Power [kJ/s = kW] 
Cp_sw = 3.9;                                                                 % Specific heat capacity[kJ/kgK], Dt_sw=10 deg C
Dt_sw_cond = 10;                                                                % Temperature difference betwen inlet and outlet of seawater cooling medium [K or Degrees C]
t_sw_cond_i = 20;
t_sw_cond_o = t_sw_cond_i + Dt_sw_cond ;
mfr_cond_sw = Q_cond/(Cp_sw*Dt_sw_cond);                                        % Mass flow rate of sea water pump [kg/s]
vfr_cond_sw = mfr_cond_sw/1025;                                              % Volume flowrate of sea water pump [m^3 /s]
eff_sw_pump =  pchip(eff_max_pump_data(:,1),eff_max_pump_data(:,2), vfr_cond_sw); % Estimation of efficiency of refrigerant pump [-]
Dp_sw_cond = 1.5;                                                            % Pressure increase of condenser sea water pump [bar]
P_pump_sw = vfr_cond_sw*Dp_sw_cond*(10^2)/eff_sw_pump;                       % Power of Sea Water pump (kW) // W = vfr*Dp
P_net = P_el - P_pump_r - P_pump_sw - P_to_pump ;

%OUTPUT
%HE:
Q_he;
UA_he;
NTU_he;
t_g_o;
t_to_he_o;

P_to_pump;
mfr_to;
vfr_to_pump;
mfr_to_ev;
mfr_to_b

%EV
Q_ev;
UA_ev;

P_pump_r;
mfr_r;
vfr_r_pump;

t_to_ev_o;
t_r_ev_o;
Dt_ev;

%TURBINE
P_el
P_net

%Condenser
Q_cond;
P_pump_sw;
mfr_cond_sw;
vfr_cond_sw;

%DRAWING
        % ORC Graph Saturation line graph
out_R245fa = Saturation_lines('R245fa');
%[x,y] = size(out_R245fa)
        % ORC Graph
t_r_cond_o; s_r_cond_o;
t_r_pump_o; s_r_pump_o; 
%s_r_pump_o_is
t_sat_p_h = CoolProp.PropsSI('T' ,'P', p_r_UP , 'Q' , 0, 'R245fa') - 273.15;
s_r_sat_l_p_h = CoolProp.PropsSI('S' ,'P', p_r_UP , 'Q' , 0, 'R245fa');
t_sat_p_h = CoolProp.PropsSI('T' ,'P', p_r_UP , 'Q' , 1, 'R245fa') - 273.15;
s_r_sat_v_p_h = CoolProp.PropsSI('S' ,'P', p_r_UP , 'Q' , 1, 'R245fa');
t_r_turb_i; s_r_turb_i;
t_r_turb_o; s_r_turb_o; 
%s_r_turb_o_is
t_sat_p_l = CoolProp.PropsSI('T' ,'P', p_r_LP, 'Q' , 1, 'R245fa') - 273.15;
s_r_sat_p_l = CoolProp.PropsSI('S' ,'P', p_r_LP , 'Q' , 1, 'R245fa');

s_r_graph = [t_r_cond_o, s_r_cond_o/(10^3) ;t_r_pump_o , s_r_pump_o/(10^3) ;t_sat_p_h, s_r_sat_l_p_h/(10^3) ;t_sat_p_h ,s_r_sat_v_p_h/(10^3) ; t_r_turb_i , s_r_turb_i/(10^3) ; t_r_turb_o, s_r_turb_o/(10^3) ; t_sat_p_l , s_r_sat_p_l/(10^3);t_r_cond_o, s_r_cond_o/(10^3);];
        % Sea water Graph
t_sw_cond_i; s_r_cond_o;
t_sw_cond_o; s_r_turb_o;
s_sw_graph = [t_sw_cond_i, s_r_cond_o/(10^3); t_sw_cond_o, s_r_turb_o/(10^3)];
        % Thermal oil Graph
t_to_ev_o; s_r_cond_o;
t_to_ev_i; s_r_turb_o;
s_to_graph = [t_to_ev_o, s_r_cond_o/(10^3);t_to_ev_i, s_r_turb_o/(10^3)];

% T/Q Graph for ORC
figure()
plot(out_R245fa(1,:),out_R245fa(2,:),'-',s_r_graph(:,2),s_r_graph(:,1),'--',s_sw_graph(:,2),s_sw_graph(:,1),'--',s_to_graph(:,2),s_to_graph(:,1),'--')
grid on
title('T/Q Graph for ORC')                                    % Title
xlabel('Percentage of heat transmission')                                                    % x-axis label
set(gca,'XTick',[])
ylabel('t_s_a_t [Degrees C]')                                                 % y-axis label
legend('R245fa Saturation line', 'ORC' ,'Sea Water', 'Themral Oil','Location','northwest')                                       % Legend

% T/Q Graph for ORC
        % Exhasut gas
t_g_i; s_r_cond_o;
t_g_o; s_r_turb_o;
s_eg_graph = [t_g_i(3), s_r_cond_o/(10^3); t_g_o, s_r_turb_o/(10^3)];
        % Thermal oil Graph
t_to_he_o; s_r_cond_o;
t_to_he_i; s_r_turb_o;
s_to_graph = [t_to_he_o, s_r_cond_o/(10^3);t_to_he_i, s_r_turb_o/(10^3)];

figure()
plot(s_eg_graph(:,2),s_eg_graph(:,1),'--',s_to_graph(:,2),s_to_graph(:,1),'--')
grid on
title('T/Q Graph for Exhaust gas heat exchanger')                                    % Title
xlabel('Percentage of heat transmission')                                                    % x-axis label
set(gca,'XTick',[])
ylabel('t_s_a_t [Degrees C]')                                                 % y-axis label
legend('Exhasut gas', 'Themral Oil','Location','northeast')        
