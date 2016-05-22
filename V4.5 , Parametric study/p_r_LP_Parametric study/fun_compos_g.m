function [ k ] = fun_compos( c, h, s, lamda)

% The composition of exhaust gas is calculated given the fuel weight
% composition in C,H,S and equivalence air/ fuel ratio (lamda)

nst_O2=(c/12.01)+((h/2.016)/2)+(s/32.066);  % moles of O2 for stoichiometric combustion
n_CO2 = c/12.01;
n_H2O = h/2.016;
n_SO2 = s/32.66;
n_O2 = (lamda-1)*nst_O2;
n_N2 = (lamda*nst_O2)*(79/21);
n_g = n_CO2+n_H2O+n_SO2+n_O2+n_N2;
m_g= ((lamda*nst_O2/0.21)*28.96)+1;
MB_g = m_g/n_g; % exhaust gas molecular weight kg/kmole

k_CO2 = ((n_CO2 /n_g) *44.01)/MB_g ;    % kg CO2/kg exhaust gas
k_H2O  = ((n_H2O /n_g) *18.016)/MB_g ;   % kg H2O/kg exhaust gas
k_SO2  = ((n_SO2 /n_g) *64.066)/MB_g ;   % kg SO2/kg exhaust gas
k_O2  = ((n_O2 /n_g) *32)/MB_g ;         % kg O2/kg exhaust gas
k_N2  = ((n_N2 /n_g )*28.016)/MB_g ;   % kg N2/kg exhaust gas

k = [k_CO2, k_H2O, k_SO2, k_O2, k_N2]; 

end

