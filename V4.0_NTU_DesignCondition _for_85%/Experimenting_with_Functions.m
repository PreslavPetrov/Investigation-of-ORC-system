eff_cor_factors_data
%Input
p_r_UP_bar =  20
p_r_LP_bar =  2
L_tg_rated = 500
El_load_est = 80
t_r_ev_o =   132.0792


%Eff_tg_rated
l = 0;
for RL = 250:250:2000
l =l+1 ;     
eff_tg_data_extr_45 = griddata(eff_tg_data(:,1),eff_tg_data(:,2),eff_tg_data(:,3), 4.5, RL, 'cubic');
eff_tg_data_extr_75 = griddata(eff_tg_data(:,1),eff_tg_data(:,2),eff_tg_data(:,3), 7.5, RL, 'cubic');
eff_tg_data_extr_150 = griddata(eff_tg_data(:,1),eff_tg_data(:,2),eff_tg_data(:,3), 15.0, RL, 'cubic');
eff_tg_data_extr (l,:) = [eff_tg_data_extr_45,eff_tg_data_extr_75,eff_tg_data_extr_150];  
end

p = [4.5 , 7.5 ,15];
n = 1;
for d = 1:1:l
Poly_tg_data_extr(d,:)=polyfit(p,eff_tg_data_extr(d,:),n);
end

ss = L_tg_rated/250;
eff_tg_rated = 0;

for i =1:1:n+1
eff_tg_rated = eff_tg_rated + Poly_tg_data_extr(ss,i)*(p_r_UP_bar)^(n-i+1);
end

%f_lcf
f_lcr = griddata(f_Lcor_st_data(:,1),f_Lcor_st_data(:,2),f_Lcor_st_data(:,3), L_tg_rated, El_load_est, 'linear');

%f_tcf
l = 0;
for t =150:25:400
l = l+1;
f_tcor_st_data_extr_45(l) = griddata(f_tcor_st_data(:,1),f_tcor_st_data(:,2),f_tcor_st_data(:,3),4.5, t, 'cubic');
f_tcor_st_data_extr_75(l) = griddata(f_tcor_st_data(:,1),f_tcor_st_data(:,2),f_tcor_st_data(:,3),7.5, t, 'cubic');
f_tcor_st_data_extr_150(l) = griddata(f_tcor_st_data(:,1),f_tcor_st_data(:,2),f_tcor_st_data(:,3), 15, t, 'cubic');
t_extr(l) =  t;
end

n = 1  ;
Poly_f_tcf_extr_45 = polyfit(t_extr,f_tcor_st_data_extr_45,n);
Poly_f_tcf_extr_75 = polyfit(t_extr,f_tcor_st_data_extr_75,n);
Poly_f_tcf_extr_150 = polyfit(t_extr,f_tcor_st_data_extr_150,n);
f_tcor_st_data_new_45 = 0;
f_tcor_st_data_new_75 = 0;
f_tcor_st_data_new_150 = 0;
for i =1:1:n+1
f_tcor_st_data_new_45 = f_tcor_st_data_new_45 + Poly_f_tcf_extr_45(i)*(t_r_ev_o)^(n-i+1);
f_tcor_st_data_new_75 = f_tcor_st_data_new_75 + Poly_f_tcf_extr_75(i)*(t_r_ev_o)^(n-i+1);
f_tcor_st_data_new_150 = f_tcor_st_data_new_150 + Poly_f_tcf_extr_150(i)*(t_r_ev_o)^(n-i+1);
end
f_tcor_st_data_pressure = [f_tcor_st_data_new_45 , f_tcor_st_data_new_75 , f_tcor_st_data_new_150 ];

Poly_f_tcf = polyfit(p,f_tcor_st_data_pressure,n);
f_tcf = 0;
for i =1:1:n+1
f_tcf = f_tcf + Poly_f_tcf(i)*(p_r_UP_bar)^(n-i+1);
end

%f_bpcf
for i =1:1:3
f_bpcf_extr_0025(i) = griddata(f_pcor_st_data(:,1),f_pcor_st_data(:,2),f_pcor_st_data(:,3),p(i), 0.025, 'cubic');
f_bpcf_extr_0050(i) = griddata(f_pcor_st_data(:,1),f_pcor_st_data(:,2),f_pcor_st_data(:,3),p(i), 0.05, 'cubic');
f_bpcf_extr_0075(i) = griddata(f_pcor_st_data(:,1),f_pcor_st_data(:,2),f_pcor_st_data(:,3),p(i), 0.075, 'cubic');
end
f_bpcf_extr_0025;
f_bpcf_extr_0050;
f_bpcf_extr_0075;

n = 1;
Poly_f_bpcf_extr_0025 = polyfit(p,f_bpcf_extr_0025,n);
Poly_f_bpcf_extr_0050 = polyfit(p,f_bpcf_extr_0050,n);
Poly_f_bpcf_extr_0075 = polyfit(p,f_bpcf_extr_0075,n);

f_bpcf_data_0025 = 0;
f_bpcf_data_0050 = 0;
f_bpcf_data_0075 = 0;
for i =1:1:n+1
f_bpcf_data_0025 = f_bpcf_data_0025 + Poly_f_bpcf_extr_0025(i)*(p_r_UP_bar)^(n-i+1);
f_bpcf_data_0050 = f_bpcf_data_0050 + Poly_f_bpcf_extr_0050(i)*(p_r_UP_bar)^(n-i+1);
f_bpcf_data_0075 = f_bpcf_data_0075 + Poly_f_bpcf_extr_0075(i)*(p_r_UP_bar)^(n-i+1);
end
f_bpcf_data = [f_bpcf_data_0025, f_bpcf_data_0050, f_bpcf_data_0075];

n = 1;
bp = [0.025,0.05,0.075];
Poly_f_bpcf = polyfit(bp,f_bpcf_data,n);
if(p_r_LP_bar >= 0.075) 
  f_bpcf = 1;
else
    f_bpcf  = 0;
    for i =1:1:n+1
     f_bpcf = f_bpcf + Poly_f_bpcf(i)*(p_r_LP_bar)^(n-i+1);
    end   
end

eff_tg_rated
f_lcr
f_tcf
f_bpcf

Eff_tg = eff_tg_rated*f_lcr*f_tcf*f_bpcf 