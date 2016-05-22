function [ cpm ] = fun_cpm_air(t1, t2)

% calculation of mean specific heat at constant pressure of air
% between t1 and t2 in degC
% t1 deg C
% t2 degC


% t1 ~= t2 
% for t1=t2 the specific heat at constant pressure in t1=t2 is calculated

B = [ 1.00186E+00 7.22153E-05 3.48066E-07 -3.21086E-10 9.30682E-14 ];
% b0-b4 for
% line 1: air

if t2 ~= t1
    cpm = 0;
    for i=1:5
     cpm = cpm + B(i)/i*(t2^i-t1^i)/(t2-t1);
    end
else % t2==t2
    cpm = 0;
    for i=1:5
     cpm = cpm + B(i)*t2^(i-1);
    end 
    
end


end

% end of function