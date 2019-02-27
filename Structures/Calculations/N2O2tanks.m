% Julia Tepper, ENAE484
% Design N2, O2 tanks for crew systems

% dimensions
    % O2 tank
    V_O2 = 110.19; %L
    V_O2 = 110.19*10^-3; % m^3
    syms r_O2;
    solve (V_O2 == 4/3*pi()*r_O2^3,r_O2)
        % r_O2 = 0.297407 m
    
    % N2 tank
    V_N2 = 73.17 * 10^-3; %m^3
    pressure = 3000; %psi
    syms r_N2;
    solve (V_N2 == 4/3*pi()*r_N2^3,r_N2)
        % r_N2 = 0.2594666 m
    
p = 2.0684e+7; % pascals
SF = 3;

% material (assume aluminum 2024)
TS = 324 * 10^6; % pascals
UTS = 469 * 10^6; % pascals

sigma = TS*SF;

t_O2 = p*r_O2/sigma %m;
t_N2 = p*r_N2/sigma %m

    % r_O2 = 0.297407 m
    % r_N2 = 0.259466 m
    % t_O2 = 0.006328772 m
    % t_N2 = 0.005521406 m