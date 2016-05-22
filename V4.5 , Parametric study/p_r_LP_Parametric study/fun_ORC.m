
function [out] = fun_ORC(t_to_he_o,p_r_LP_bar,p_r_UP_bar,vfr_r,UA_ev_d,mfr_to,t_to_he_i,Refrigerant,vfr_sw,P_to_pump,L_tg_rated)
thermal_properties_data_to
eff_pump_data
Gen_mech_eff_data

%2. Interm Loop                                                
t_to_ev_i = t_to_he_o; % Inlet Temperature of Thermal oil to Evaporator [Degrees C]

%4. ORC
%PUMP                                                           
p_r_LP = p_r_LP_bar*100000;                                                  % Lower Pressure of ORC system [Pa], needed for CoolProp
p_r_UP = p_r_UP_bar*100000;                                                  % Upper Pressure of ORC system [Pa], needed for CoolProp

tK_r_cond_sat = CoolProp.PropsSI('T' ,'P', p_r_LP, 'Q' , 0, 'R245FA');       % Saturation temeprature of Refrigerant at Lower pressure, p_r_LP [K]
t_r_cond_sat = tK_r_cond_sat - 273.15;                                       % Saturation temeprature of Refrigerant at Lower pressure, p_r_LP [Degrees C]
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

mfr_r = vfr_r*dens_r_pump ;                                                 % Requred Mass flowrate of Refrigerant [m^3/s]
eff_r_pump =  pchip(eff_max_pump_data(:,1),eff_max_pump_data(:,2), vfr_r);  % Estimation of efficiency of refrigerant pump [-]
P_pump_r = vfr_r*(p_r_UP-p_r_LP)/(eff_r_pump*(10^3));  

%EVAPORATOR
eff_ev_est = 0.8;                %FIRST estiamte
t_to_ev_o_est = 120;             %FIRST estiamte
t_r_ev_i = t_r_pump_o;
h_r_ev_i = h_r_pump_o/(10^3);   %in kJ/kg
n = 20;
i = 1;
while i < n
Cp_to_ev = pchip(TO_Properties(:,1),TO_Properties(:,2), (t_to_ev_o_est+t_to_ev_i)/2);
Cp_to_ev_min = pchip(TO_Properties(:,1),TO_Properties(:,2), (t_r_ev_i+t_to_ev_i)/2);
t_to_ev_o = t_to_ev_i - (Cp_to_ev_min/Cp_to_ev)*(t_to_ev_i-t_r_ev_i)*eff_ev_est;
t_to_ev_o_est = t_to_ev_o;

Cp_to_he_i = pchip(TO_Properties(:,1),TO_Properties(:,2), t_to_he_i);
Cp_to_he_o = pchip(TO_Properties(:,1),TO_Properties(:,2), t_to_he_o);
Cp_to_ev_o = pchip(TO_Properties(:,1),TO_Properties(:,2), t_to_ev_o);
mfr_to_ev = mfr_to*(Cp_to_he_i*t_to_he_i - Cp_to_he_o*t_to_he_o)/(Cp_to_ev_o*t_to_ev_o - Cp_to_he_o*t_to_he_o);
    if(mfr_to_ev > mfr_to)
        mfr_to_ev = mfr_to;
    end
mfr_to_b = mfr_to - mfr_to_ev;

%mfr_to*Cp_to_he_i*t_to_he_i - mfr_to_b*Cp_to_he_o*t_to_he_o -mfr_to_ev*Cp_to_ev_o*t_to_ev_o 

h_r_ev_o = h_r_ev_i + mfr_to_ev*Cp_to_ev*(t_to_ev_i - t_to_ev_o)/mfr_r;  %kJ/kg

t_r_ev_o = CoolProp.PropsSI('T' ,'P', p_r_UP, 'H' , h_r_ev_o*(10^3), 'R245FA') -273.15;
h_r_ev_sat = CoolProp.PropsSI('H' ,'P', p_r_UP, 'Q' , 1, 'R245FA')/(10^3);
t_r_ev_sat = CoolProp.PropsSI('T' ,'P', p_r_UP, 'Q' , 1, 'R245FA') -273.15;
Dt_ev = t_r_ev_o - t_r_ev_sat;

% Correction for UA for off-design operation 
% mfr_to_ev_d  = 10.1565 kg/s
UA_ev_p = UA_ev_d*(mfr_to_ev/10.1565)^0.65;

NTU_ev = UA_ev_p/Cp_to_ev*mfr_to_ev;
eff_ev = 1-exp(-NTU_ev);
eff_ev_est = eff_ev;

    if h_r_ev_o < h_r_ev_sat 
     t_to_he_i = t_to_he_i -1;
     n = n+1;
    end
i = i + 1;
end

Q_ev = mfr_to_ev*Cp_to_ev*(t_to_ev_i - t_to_ev_o);

%TURBINE
h_r_turb_i = h_r_ev_o;                                                       % Inlet Enthalpy of Refrigerant at Turbine inlet [J/kg]
t_r_turb_i = t_r_ev_o;
tK_r_turb_i = t_r_turb_i +273.15;
s_r_turb_i = CoolProp.PropsSI('S' ,'P', p_r_UP, 'T' , tK_r_turb_i, Refrigerant);  % Inlet Entropy of Refrigerant at Turbine [J/kgK]                                                  % 
s_r_turb_o_is = s_r_turb_i;                                                  % Isentropic Outlet Entropy of Refrigerant at Turbine [J/kgK]    
h_r_turb_o_is = CoolProp.PropsSI('H' ,'P', p_r_LP, 'S' , s_r_turb_o_is, Refrigerant);

% Turb Isentropic efficiency
N_turb = 20000 ; % RPM of turbine in rev per min
Ns_turb = N_turb/60 ; % RPS of turbine in rev per sec
Dens_r_turb_o_is = CoolProp.PropsSI('D' ,'P', p_r_LP, 'H' , h_r_turb_o_is, Refrigerant);
vfr_r_turb_o = mfr_r/(Dens_r_turb_o_is);
delh_is = (h_r_turb_i*(10^3)) - h_r_turb_o_is;
Ns = Ns_turb*((vfr_r_turb_o)^0.5)/(delh_is^(0.75));
eff_turb_is = fun_is_eff_turb(Ns) ;                                                           % Isentropic Efficiency of Turbine

h_r_turb_o = h_r_turb_i - eff_turb_is*(h_r_turb_i - h_r_turb_o_is/(10^3));          % Outlet Entropy of Refrigerant at Turbine [J/kgK]   
t_r_turb_o = CoolProp.PropsSI('T' ,'P', p_r_LP, 'H' , h_r_turb_o*(10^3), Refrigerant) -273.15;
s_r_turb_o = CoolProp.PropsSI('S' ,'P', p_r_LP, 'H' , h_r_turb_o*(10^3), Refrigerant);

El_load_est = 80;
for it =1:1:5
eff_m  = fun_mech_eff_turb (L_tg_rated,El_load_est);
P_el = mfr_r*(h_r_turb_i - h_r_turb_o)*eff_m ;                            % Isentropic Power of Turbine (REfrigerant power) [kJ/s = kW]
El_load_est =100*P_el/L_tg_rated;
end

Dens_r_turb_o = CoolProp.PropsSI('D' ,'P', p_r_LP, 'H' , h_r_turb_o*(10^3), Refrigerant);
vfr_r_turb_o = mfr_r/(Dens_r_turb_o);
VH = ((vfr_r_turb_o)^(0.5))/((delh_is)^(0.25));   % Size parameter [m]

Dens_r_turb_o = CoolProp.PropsSI('D' ,'P', p_r_UP, 'H' , h_r_turb_i*(10^3), Refrigerant);
vfr_r_turb_i = mfr_r/Dens_r_turb_o;
VR = vfr_r_turb_o/vfr_r_turb_i;   % Voluem exp ratio

%CONDENSER
Q_cond =  mfr_r*(h_r_turb_o - h_r_cond_o/(10^3)) ;                           % Condenser Thermodynamic Power [kJ/s = kW] 
Cp_sw = 3.9;                                                                 % Specific heat capacity[kJ/kgK], Dt_sw=10 deg C    
t_sw_cond_i = 20;
mfr_sw= vfr_sw*1025;
t_sw_cond_o = t_sw_cond_i+ Q_cond/(mfr_sw*(Cp_sw));                                        % Mass flow rate of sea water pump [kg/s]                                              % Volume flowrate of sea water pump [m^3 /s]
eff_sw_pump =  pchip(eff_max_pump_data(:,1),eff_max_pump_data(:,2), vfr_sw); % Estimation of efficiency of refrigerant pump [-]
Dp_sw_cond = 1.5;                                                            % Pressure increase of condenser sea water pump [bar]
P_pump_sw = vfr_sw*Dp_sw_cond*(10^2)/eff_sw_pump;                       % Power of Sea Water pump (kW) // W = vfr*Dp
P_net = P_el - P_pump_r - P_pump_sw - P_to_pump ;
P_pumps = P_pump_r + P_pump_sw + P_to_pump;


%OUTPUT
out(1,:) = P_pump_r;
out(2,:) = h_r_ev_o;
out(3,:) = h_r_ev_sat;
out(4,:) = t_to_he_i;
out(5,:) = Q_ev;
out(6,:) = eff_turb_is;
out(7,:) = P_pumps;
out(8,:) = t_sw_cond_o;
out(9,:) = Q_cond;
out(10,:) = mfr_to_b; 
out(11,:) = mfr_to_ev;
out(12,:) = P_el ;
out(13,:) = P_net ;
out(14,:) = eff_m;

end