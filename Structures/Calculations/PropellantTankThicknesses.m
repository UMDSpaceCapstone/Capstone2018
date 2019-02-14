% Julia Tepper, ENAE848
% Look at thicknesses for propellant tanks

% dimensions
    % N204 tank stage 2
    rN204_2 = 0.415; % m
    cN204_2 = 0.208; % m
    hN204_2 = 2.668; % m
    nN204_2 = 2;
    vN204_2 = (4/3)*pi()*rN204_2^2*cN204_2+pi()*rN204_2^2*hN204_2; % m^3
    
    % MMH tank stage 2
    rMMH_2 = 0.415; % m
    cMMH_2 = 0.208; % m
    hMMH_2 = 2.276; % m
    nMMH_2 = 2; % m
    vMMH_2 = (4/3)*pi()*rMMH_2^2*cMMH_2+pi()*rMMH_2^2*hMMH_2; % m^3
    
    % N204 tank stage 1
    rN204_1 = 0.828; % m
    cN204_1 = 0.414; % m
    hN204_1 = 1.834; % m
    nN204_1 = 2;
    vN204_1 = ((4/3)*pi()*rN204_1^2*cN204_1)+(pi()*rN204_1^2*hN204_1); % m^3

    % MMH tank stage 1
    rMMH_1 = 0.777; % m
    cMMH_1 = 0.389; % m
    hMMH_1 = 1.832; % m
    nMMH_1 = 2; % m
    vMMH_1 = (4/3)*pi()*rMMH_1^2*cMMH_1+pi()*rMMH_1^2*hMMH_1; % m^3

% pressure (for all tanks)
p = 96526.6; % pascals (N/m^2)
 
% material (assume aluminum 2024)
TS = 324 * 10^6; % pascals
UTS = 469 * 10^6; % pascals

SF = 3; % safety factor

sigma = TS*SF;

% N204 2
    
    tN204_2 = p*rN204_2/sigma % m
    
% MMH 2

    tMMH_2 = p*rMMH_2/sigma % m
    
% N204 1

    tN204_1 = p*rN204_1/sigma % m
    
% MMH 1

    tMMH_1 = p*rMMH_1/sigma % m
    



