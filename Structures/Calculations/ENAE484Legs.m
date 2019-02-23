% Juliette Abbonizio and Mo Khaial Lander Legs Calculations
% using the mass and the delta vs calculate the forces
% assuming the legs will be made out of 6061 aluminum

clear all
% masses (kg)
mpr_des = .10*27312;
mi_des = 2242.4; %includes structural mass
mtot_asc = 13264;
mtot = mpr_des + mtot_asc + mi_des-490;

%moon gravity (m/s^2)
g = 1.62;

% velocity at touchdown (m/s)
v = 2;
% time after touchdown (s)
t = 1;
%impulse (N-s)
%I = mtot*v;

%weight
W = mtot*g;

%energy
Ke = 1/2*mtot*v^2;

% yield strength (sigma) (Pa)
%O = Pcr/A
% O = 276e6;
% P = O.*A;

% vertical acceleration loading
%diameter and thickness for a circular ring
 do = .05:.01:.15;
 dt = .01:.001:.041;
 count = 1;
 count1 = 1;
for i = 1:length(do)
    for j = 1:length(dt)
        di(i,j) = do(i)-dt(j);
        % length of the legs
        L = 1.5:.1:1.8;

        % use area of an annulus (ring)
        A(i,j) = pi/4.*(do(i).^2-di(i,j).^2); % circular ring

        %moment of intertia of a ring (m^4)
        I(i,j) = pi/64.*(do(i).^4-di(i,j).^4);

        %E for 6061 aluminum (Pa)
        E_al = 68.9e9;
        
        %Density of 6061 aluminum (kg/m^3)
        rho_al = 2710;
        
        % ultimate safety factor
        SFu = 2;

        %buckling (using aluminum)
        %Pcr = (pi^2*E*I)/(KL)^2
        Pcr = (mtot*(v/t)*SFu); %max buckling kN
        K = .5; % (fixed-fixed)
        for k = 1:length(L)
            Pcrb_al(i,j) = ((pi^2*E_al.*I(i,j))./(K*L(k))^2); %actual buckling kN (Aluminum)
            
            % Volume of Hollow Cylinder V = pi/4*h*(do^2-di^2)
            V = pi/4*L(k)*(do(i)^2-(do(i)-dt(j))^2);
            
            % Mass of 1 leg (kg)
            M_leg_al = rho_al*V;
             
            if Pcrb_al(i,j) > Pcr
                Pact_al(1,count) = Pcrb_al(i,j);
                Pact_al(2,count) = do(i);
                Pact_al(3,count) = dt(j);
                Pact_al(4,count) = L(k);
                Pact_al(5,count) = 4*M_leg_al;
                count = count+1;
            end
       
        end
    end
    
  
end

%upper leg
%force; outer diameter(m) ; thickness(m) ; length(m) ; weight (kg)
[1760093.9446; 0.1; 0.013; 1.8; 37.25]; %first try
[2433045.305; 0.1; 0.02; 1.8; 55.17] %second try

%lower leg

% horizontal acceleration loading

