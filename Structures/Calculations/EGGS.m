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
ri_ch4 = (V_ch4/((4/3)*pi))^(1/3);
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

MS_ch4 = Sy/43238596-1;
% stress, thickness, mass
% [43238596, 0.003, 84]

% LOX tank
P_lox = 206.8e3;
V_lox = 10.48;
ri_lox = (V_lox/((4/3)*pi))^(1/3);
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

MS_lox = Sy/40106064-1;
% stress, thickness, mass
% [40106064, 0.0035, 114]

t2 = .005:.001:.01;

% MMH tank
P_mmh = 10e6; %Pa
V_mmh = .07 + .1798; %m3
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

MS_mmh = Sy/390692080-1;
% stress, thickness, mass
% [390692080, 0.005, 27]

% N2O4 tank
P_no = 10e6;
V_no = .092 + .2373;
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

MS_no = Sy/428385116-1;

% stress, thickness, mass
% [428385116, 0.005, 32.4]