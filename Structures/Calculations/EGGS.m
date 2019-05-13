% Juliette Abbonizio EGGS
%% total mass
mpr = 1065; %kg
mtot = 18999.20; %kg

%% tank calcs
% CH4 tank
P_ch4 = 206.8e3;
V_ch4 = 8.27;
%V = 4/3*pi*r^3
t1 = .001:.0005:.004;
Sy = 503e6; %MPa %yield for aluminum
count = 1;
ri_ch4 = (V_ch4/((2/3)*pi))^(1/3);
rho_al = 2810;
SA_ch4 = 2*pi*ri_ch4^2;
for i = 1:length(t1)
    sig_ch4(i) = (P_ch4*ri_ch4)/(2*t1(i)); %axial stress;
    m_ch4(i) = rho_al*SA_ch4*t1(i);
    if sig_ch4(i) < Sy
        Sig_act_ch4(count,1) = sig_ch4(i);
        Sig_act_ch4(count,2) = t1(i);
        Sig_act_ch4(count,3) = m_ch4(i);
        count = count+1;
    end
end
ro_ch4 = .001+ri_ch4;
% stress, thickness, mass
% [54477217, 0.003, 132.3]

% LOX tank
P_lox = 206.8e3;
V_lox = 10.48;
ri_lox = (V_lox/((2/3)*pi))^(1/3);
rho_al = 2810;
SA_lox = 2*pi*ri_lox^2;
count2 = 1;
for i = 1:length(t1)
    sig_lox(i) = (P_lox*ri_lox)/(2*t1(i)); %axial stress;
    m_lox(i) = rho_al*SA_lox*t1(i);
    if sig_lox(i) < Sy
        Sig_act_lox(count2,1) = sig_lox(i);
        Sig_act_lox(count2,2) = t1(i);
        Sig_act_lox(count2,3) = m_lox(i);
        count2 = count2 +1;
    end
end
ro_lox = ri_lox + .001;
% stress, thickness, mass
% [5053047, 0.0035, 181]

t2 = .005:.001:.01;

% MMH tank
P_mmh = 10e6; %Pa
V_mmh = .07; %m3
ri_mmh = (V_mmh/((4/3)*pi))^(1/3);
rho_al = 2810;
SA_mmh = 4*pi*ri_mmh^2;
count2 = 1;
for i = 1:length(t2)
    sig_mmh(i) = (P_mmh*ri_mmh)/(2*t2(i)); %axial stress;
    m_mmh(i) = rho_al*SA_mmh*t2(i);
    if sig_mmh(i) < Sy
        Sig_act_mmh(count2,1) = sig_mmh(i);
        Sig_act_mmh(count2,2) = t2(i);
        Sig_act_mmh(count2,3) = m_mmh(i);
        count2 = count2 +1;
    end
end
ro_mmh = ri_mmh + .001;

% stress, thickness, mass
% [255664135, 0.005, 11.5]

% N2O4 tank
P_no = 10e6;
V_no = .092;
ri_no = (V_no/((4/3)*pi))^(1/3);
rho_al = 2810;
SA_no = 4*pi*ri_no^2;
count2 = 1;
for i = 1:length(t2)
    sig_no(i) = (P_no*ri_no)/(2*t2(i)); %axial stress;
    m_no(i) = rho_al*SA_no*t2(i);
    if sig_no(i) < Sy
        Sig_act_no(count2,1) = sig_no(i);
        Sig_act_no(count2,2) = t2(i);
        Sig_act_no(count2,3) = m_no(i);
        count2 = count2 +1;
    end
end
ro_no = ri_no + .001;

% stress, thickness, mass
% [280048385, 0.005, 13.8]