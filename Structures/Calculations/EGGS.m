% Juliette Abbonizio EGGS
%% total mass
mpr = 1034.20; %kg
mtot = 18474; %kg

%% tank calcs
% CH4 tank
P_ch4 = 206.8e3;
V_ch4 = 8.27;
%V = 4/3*pi*r^3
t = .001:.0005:.004;
Sy = 503e6; %MPa %yield for aluminum
count = 1;
ri_ch4 = (V_ch4/((2/3)*pi))^(1/3);
rho_al = 2810;
SA_ch4 = 2*pi*ri_ch4^2;
for i = 1:length(t)
    sig(i) = (P_ch4*ri_ch4)/(2*t(i)); %axial stress;
    m_ch4(i) = rho_al*SA_ch4*t(i);
    if sig(i) < Sy
        Sig_act_ch4(count,1) = sig(i);
        Sig_act_ch4(count,2) = t(i);
        Sig_act_ch4(count,3) = m_ch4(i);
        count = count+1;
    end
end
ro_ch4 = .001+ri_ch4;
% stress, thickness, mass
% [40857913, 0.004, 176.4]
% [163431652, 0.001, 44]

% LOX tank
P_lox = 206.8e3;
V_lox = 10.48;
ri_lox = (V_lox/((2/3)*pi))^(1/3);
rho_al = 2810;
SA_lox = 2*pi*ri_lox^2;
count2 = 1;
for i = 1:length(t)
    sig(i) = (P_lox*ri_lox)/(2*t(i)); %axial stress;
    m_lox(i) = rho_al*SA_lox*t(i);
    if sig(i) < Sy
        Sig_act_lox(count2,1) = sig(i);
        Sig_act_lox(count2,2) = t(i);
        Sig_act_lox(count2,3) = m_lox(i);
        count2 = count2 +1;
    end
end
ro_lox = ri_lox + .001;
% stress, thickness, mass
% [176856664, 0.001, 52]

%% skirt vertical force (buckling)
%solve for buckling and stress in the shear panels (use mass of whole
%vehicle)
F1 = 6*9.81*mtot;

%% skirt horizontal force (bending and shear)
F2 = 2*9.81*mtot;