% Juliette Abbonizio Lander Legs Calculations
% using the mass and the delta vs calculate the forces
% assuming the legs will be made out of 6061 aluminum

clear all
% masses (kg)
mpr_des = .10*22745;
mi_des = 2242.4; %includes structural mass
mtot_asc = 13015;
mtot = mpr_des + mtot_asc + mi_des;

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
 do = .05:.01:.5;
 dt = .001:.001:.041;
for i = 1:length(do)
    for j = 1:length(dt)
        di(i,j) = do(i)-dt(j);
        % length of the legs
        L = 1.8;

        % use area of an annulus (ring)
        A(i,j) = pi/4.*(do(i).^2-di(i,j).^2); % circular ring

        %moment of intertia of a ring (m^4)
        I(i,j) = pi/64.*(do(i).^4-di(i,j).^4);

        %E for 6061 aluminum (Pa)
        E = 68.9e9;

        % ultimate safety factor
        SFu = 2;

        %buckling (using aluminum)
        %Pcr = (pi^2*E*I)/(KL)^2
        Pcr = (mtot*(v/t)*SFu); %max buckling kN
        K = .5; % (fixed-fixed)
        Pcrb(i,j) = ((pi^2*E.*I(i,j))./(K*L)^2); %actual buckling kN
        
        
    end
    
  
end

count = 1;
for i = 1:length(do)
    for j = 1:length(dt)
        if Pcrb(i,j) < Pcr
            Pact(1,count) = Pcrb(i,j);
            Pact(2,count) = i;
            Pact(3,count) = j;
            count = count+1;
        end
    end
end


% horizontal acceleration loading

