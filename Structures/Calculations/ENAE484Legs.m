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

% vertical velocity at touchdown (m/s)
vv = 2;
% time after touchdown (s)
t = 1;
%impulse (N-s)
%I = mtot*v;

%weight
W = mtot*g;

%energy
Ke = 1/2*mtot*vv^2;

% yield strength (sigma) (Pa)
% O = Pcr/A
% O = 276e6;
% P = O.*A;

% vertical acceleration loading
% diameter and thickness for a circular ring
 do = .1:.01:.25;
 dt = .03:.001:.060;
 count = 1;
 count1 = 1;
for i = 1:length(do)
    for j = 1:length(dt)
        di(i,j) = do(i)-dt(j);
        ro(i) = do(i)/2;
        ri(i,j) = di(i,j)/2; 
        % length of the legs
        L = 1.5:.1:3;

        % use area of an annulus (ring)
        A(i,j) = pi/4.*(do(i).^2-di(i,j).^2); % circular ring

        %moment of intertia of a ring (m^4)
        I(i,j) = pi/64.*(do(i).^4-di(i,j).^4);

        %E for 6061 aluminum (Pa) %7075
        E_al = 68.9e9;
        
        %Density of 6061 aluminum (kg/m^3)
        rho_al = 2710;
        
        % ultimate safety factor
        SFu = 2;

        %FORCES
        %Pcr = (pi^2*E*I)/(KL)^2
        F = (mtot*(vv/t)*SFu); % estimated landing force
        K = .5; % (fixed-pinned)
        for k = 1:length(L)
            Pcr(i,j) = ((pi^2*E_al.*I(i,j))./(K*L(k))^2); %Buckling N
            Sigma(i,j) = F/A(i,j); %Normal Stress (Pa)
           
            % Volume of Hollow Cylinder V = pi/4*h*(do^2-di^2)
            V = pi/4*L(k)*(do(i)^2-(do(i)-dt(j))^2);
            
            % Mass of 1 leg (kg)
            M_leg_al = rho_al*V;
             
            if Pcr(i,j) > F
                Pact(1,count) = Pcr(i,j);
                Pact(2,count) = do(i);
                Pact(3,count) = dt(j);
                Pact(4,count) = L(k);
                Pact(5,count) = 4*M_leg_al;
                count = count+1;
            end
       
        end
    end
    
  
end

% horiztonal acceleration loading
vh = 
do = .1:.01:.25;
dt = .03:.001:.060;
count = 1;
count1 = 1;
for i = 1:length(do)
    for j = 1:length(dt)
        di(i,j) = do(i)-dt(j);
        ro(i) = do(i)/2;
        ri(i,j) = di(i,j)/2; 
        % length of the legs
        L = 1.5:.1:3;

        % use area of an annulus (ring)
        A(i,j) = pi/4.*(do(i).^2-di(i,j).^2); % circular ring

        %moment of intertia of a ring (m^4)
        I(i,j) = pi/64.*(do(i).^4-di(i,j).^4);
        
        %Density of 6061 aluminum (kg/m^3)
        rho_al = 2710;
        
        % ultimate safety factor
        SFu = 2;

        %FORCES
        %Pcr = (pi^2*E*I)/(KL)^2
        F = (mtot*(vh/t)*SFu); % estimated landing force
        for k = 1:length(L)
            Tau = ((4*V)/(3*A))*((ro(i,j)^2+(ro(i,j)*ri(i,j))+ri(i,j)^2)/(ro(i,j)^2+ri(i,j)^2));
            M = F*L(k);
            c = ro(i,j);
            Sigmab(i,j) = (M*c)/I(i,j); %Bending Stress (Pa)
            
            % Volume of Hollow Cylinder V = pi/4*h*(do^2-di^2)
            V = pi/4*L(k)*(do(i)^2-(do(i)-dt(j))^2);
            
            % Mass of 1 leg (kg)
            M_leg_al = rho_al*V;
             
            if Tau(i,j) > F || Sigmab(i,j) > F
                Pact(1,count) = Sigmab(i,j);
                Pact(2,count) = Tau(i,j);
                Pact(3,count) = do(i);
                Pact(4,count) = dt(j);
                Pact(5,count) = L(k);
                Pact(6,count) = 4*M_leg_al;
                count = count+1;
            end
       
        end
    end
    
  
end