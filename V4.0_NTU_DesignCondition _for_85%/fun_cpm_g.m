function [ cpm ] = fun_cpm(t1, t2, k)


% calculation of mean specific heat at constant pressure of exhaust gas
% between t1 and t2 in degC
% t1 deg C
% t2 degC
% k(1:5) weighted fraction of exhaust gas in CO2, H2O, SO2, O2, N2 kg/kg of
% exhaust gas

% t1 ~= t2 
% for t1=t2 the specific heat at constant pressure in t1=t2 is calculated

B = [   8.28204E-01	9.81404E-04	-7.90052E-07  3.28413E-10 -5.46602E-14
        1.85042E+00	2.88423E-04	 7.14063E-07 -4.78786E-10  9.43951E-14
        5.92914E-01	6.38217E-04	-6.18659E-07  2.83124E-10 -4.91597E-14  
        9.02430E-01	3.61332E-04	-1.64362E-07  2.16244E-11  3.54211E-15
        1.03693E+00	2.78472E-05	 3.92958E-07 -3.13739E-10  7.20044E-14 ];
% b0-b4 for
% line 1: CO2
% line 2: H2O
% line 3: SO2
% line 4: O2
% line 5: N2

b=zeros(1,5);
for i =1:5
    b(i) = k(1)*B(1,i)+k(2)*B(2,i)+k(3)*B(3,i)+k(4)*B(4,i)+k(5)*B(5,i);
end

if t2 ~= t1
    cpm = 0;
    for i=1:5
     cpm = cpm + (b(i)/i)*(t2^i-t1^i)/(t2-t1);
    end
else % t2==t1
    cpm = 0;
    for i=1:5
     cpm = cpm + b(i)*t2^(i-1);
    end 
    
end


end

% end of function