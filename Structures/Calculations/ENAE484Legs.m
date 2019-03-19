% Juliette Abbonizio and Mo Khaial Lander Legs Calculations
% using the mass and the delta vs calculate the forces
% assuming the legs will be made out of 7075 aluminum

clear all
% masses (kg)
mpr = 31426*.9;
mt = 40528;
mtot = mt - mpr;

%moon gravity (m/s^2)
g = 1.62;

% vertical velocity at touchdown (m/s)
vv = 2;
% time after touchdown (s)
t = 1;

%weight
W = mtot*g;

%energy
Ke = 1/2*mtot*vv^2;

% vertical acceleration loading
% diameter and thickness for a circular ring
 do = .3:.01:.5;
 dt = .005:.001:.01;
 count = 1;
 count1 = 1;
for i = 1:length(do)
    for j = 1:length(dt)
        di(i,j) = do(i)-dt(j);
        ro(i) = do(i)/2;
        ri(i,j) = di(i,j)/2; 
        % length of the legs
        L = 2:.1:3;

        % use area of an annulus (ring)
        A(i,j) = pi*(ro(i)^2-ri(i,j)^2); % circular ring

        %moment of intertia of a ring (m^4)
        I(i,j) = pi/4.*(ro(i).^4-ri(i,j).^4);

        %E for 7075 aluminum (Pa)
        E = 71.7e9;
        
        %Density of 7075 aluminum (kg/m^3)
        rho_al = 2810;
        
        % ultimate safety factor
        SFu = 2;

        %FORCES
        F = (mtot*(vv/t)*SFu); % estimated landing force
        K = .7; % (fixed-pinned)
        for k = 1:length(L)
            Pcr(i,j) = ((pi^2*E.*I(i,j))./(K*L(k))^2); %Buckling N
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

figure(1)
plot(dt,Pcr)
title('Thickness vs. Buckling')
xlabel('thickness')
ylabel('buckling')
legend('length = 2', 'length = 2.1', 'length = 2.2', 'length = 2.3','length = 2.4', 'length = 2.5', 'length = 2.6', 'length = 2.7', 'length = 2.8','length = 2.9', 'length = 3')

% horiztonal acceleration loading
vh = 1;
do = .1:.01:.25;
dt = .03:.001:.06;
count = 1;
count1 = 1;
for i = 1:length(do)
    for j = 1:length(dt)
        di(i,j) = do(i)-dt(j);
        ro(i) = do(i)/2;
        ri(i,j) = di(i,j)/2; 
        % length of the legs
        L = 1:.1:3;

        % use area of an annulus (ring)
        A(i,j) = pi/4.*(do(i).^2-di(i,j).^2); % circular ring

        %moment of intertia of a ring (m^4)
        I(i,j) = pi/64.*(do(i).^4-di(i,j).^4);
        
        %Density of 6061 aluminum (kg/m^3)
        rho_al = 2810;
        
        % ultimate safety factor
        SFu = 2;

        %FORCES
        %Pcr = (pi^2*E*I)/(KL)^2
        F = (mtot*(vh/t)*SFu); % estimated landing force
        Vf = F; % estimated shear force %is this right???
        
        for k = 1:length(L)
            Tau(i,j) = ((4*Vf)/(3*A(i,j)))*((ro(i)^2+(ro(i)*ri(i,j))+ri(i,j)^2)/(ro(i)^2+ri(i,j)^2)); %shear stress
            %Tau(i,j) = T*r/J %Torsion or shear stress?? who knows
            M = F*L(k)/2;
            c(i,j) = ro(i);
            Sigmab(i,j) = (M*c(i,j))/I(i,j); %Bending Stress (Pa)
            
            % Volume of Hollow Cylinder V = pi/4*h*(do^2-di^2)
            V = pi/4*L(k)*(do(i)^2-(do(i)-dt(j))^2);
            
            % Mass of 1 leg (kg)
            M_leg_al = rho_al*V;
             
            if Tau(i,j) > F && Sigmab(i,j) > F
                Fact(1,count) = Sigmab(i,j);
                Fact(2,count) = Tau(i,j);
                Fact(3,count) = do(i);
                Fact(4,count) = dt(j);
                Fact(5,count) = L(k);
                Fact(6,count) = 4*M_leg_al;
                count = count+1;
            end
       
        end
    end
    
  
end