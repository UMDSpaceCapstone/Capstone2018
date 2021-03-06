% Juliette Abbonizio and Mo Khaial Lander Legs Calculations
% using the mass and the delta vs calculate the forces

clear all
% masses (kg)
mpr = 23944+7482;
mpr_des = .50*mpr;
mtot = 40528-(mpr-mpr_des);

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

% vertical acceleration loading
% diameter and thickness for a circular ring
 do = .1:.01:.25;
 dt = .004:.001:.007;
 count = 1;
 count1 = 1;
for i = 1:length(do)
    for j = 1:length(dt)
        di(i,j) = do(i)-dt(j);
        ro(i) = do(i)/2;
        ri(i,j) = di(i,j)/2; 
        % length of the legs
        L = 1.8:.1:2.5;

        % use area of an annulus (ring)
        A(i,j) = pi/4*(do(i)^2-di(i,j)^2); % circular ring

        %moment of intertia of a ring (m^4)
        I(i,j) = pi/64*(do(i)^4-di(i,j)^4);

        %E for 7075 aluminum (Pa) 
        E = 71.7e9;
        
        %Density of 7075 aluminum (kg/m^3)
        rho_al = 2810;
        
        % ultimate safety factor
        SFu = 2;

        %FORCES
        %Pcr = (pi^2*E*I)/(KL)^2
        F = (mtot*(vv/t)*SFu); % estimated landing force
        K = .65; % (fixed-pinned)
        for k = 1:length(L)
            Pcr(i,j) = ((pi^2*E*I(i,j))/(K*L(k))^2); %Buckling N
            Sigma(i,j) = F/A(i,j); %Normal Stress (Pa)
           
            % Volume of Hollow Cylinder V = pi/4*h*(do^2-di^2)
            V = pi/4*L(k)*(do(i)^2-(do(i)-dt(j))^2);
            
            % Mass of 1 leg (kg)
            M_leg_al = rho_al*V;

            if Pcr(i,j) > F/2 &&  4*M_leg_al < 200 && 4*M_leg_al > 40
                Pact(count,1) = Pcr(i,j);
                Pact(count,2) = do(i);
                Pact(count,3) = dt(j);
                Pact(count,4) = L(k);
                Pact(count,5) = 4*M_leg_al;
                count = count+1;
            end
            
            
       
        end
    end
    
  
end
% [force, outer diameter, thickness, length, weight]
% [4368029, 0.2, 0.007, 2, 48.57] %update CAD
Sy = 503e6;
MS = (4368029/F)-1
%% horiztonal acceleration loading
% masses (kg)
mpr = 23944+7482;
mpr_des = .50*mpr;
mi_des = 2242.4; %includes structural mass
mtot = 40528-(mpr-mpr_des);

vh = 1; %horizontal velocity (m/s)
t = 1; %time to touchdown (s)
P = 99260/2; % vertical force
W = (.8*P)/2; %horizontal force (N) using friction (.8) coefficient of friction
x = 0:.1:1.2; %distance from the end of the beam to wherever we want to analyze
a = 0; % W occurs at the tip
l = 2;
E = 71.7e9;
do = .2;
dt = .007;
di = do-dt;
ro = do/2;
ri = di/2;
I = pi/64*(do^4-di^4);
A = pi/4*(do^2-di^2);
    
% beams under axial compression and transverse loading
% MAKE GRAPHS WITH DIFFERENT A VALUES ASSUMING MAX BENDING WHEN a = 0
% table 8.8 1a

%constants
k = sqrt(P/(E*I));
F1 = cos(k*x);
F2 = sin(k*x);
F3 = 1-cos(k*x);
F4 = k*x-sin(k*x);

%constants
C1 = cos(k*l);
C2 = sin(k*l);
C3 = 1-cos(k*l);
C4 = k*l-sin(k*l);
Ca1 = cos(k*(l-a));
Ca2 = sin(k*(l-a));
Ca3 = 1-cos(k*(l-a));
Ca4 = k*(l-a)-sin(k*(l-a));
Ca5 = k^2/2*(l-a)^2-Ca3;
Ca6 = k^3/6*(l-a)^3-Ca4;

Ra = 0;
Siga = (W/P)*(Ca3/C1);
Ma = 0;
ya = (-W/(k*P))*((C2*Ca3-C1*Ca4)/C1);
yb = 0;
Rb = W;
Sigb = 0;
Mb = -(W/k)*((C2*Ca3+C1*Ca2)/C1);
for i = 1:length(x)
    if (i/10) < a
        Fa1(i) = 0*cos(k*(x(i)-a));
        Fa2(i) = sin(k*(x(i)-a));
        Fa3(i) = 0*(1-cos(k*(x(i)-a)));
        Fa4(i) = k*(x(i)-a)-sin(k*(x(i)-a));
        Fa5(i) = ((k^2)/2)*((x(i)-a)^2)-Fa3(i);
        Fa6(i) = ((k^3)/6)*((x(i)-a)^3)-Fa4(i);
    else
        Fa1(i) = cos(k*(x(i)-a));
        Fa2(i) = sin(k*(x(i)-a));
        Fa3(i) = 1-cos(k*(x(i)-a));
        Fa4(i) = k*(x(i)-a)-sin(k*(x(i)-a));
        Fa5(i) = ((k^2)/2)*((x(i)-a)^2)-Fa3(i);
        Fa6(i) = ((k^3)/6)*((x(i)-a)^3)-Fa4(i);
    end
end
V = (Ra.*F1)-(Ma.*k.*F2)-(Siga.*P.*F1)-(W.*Fa1); %transverse shear
M = (Ma.*F1)+((Ra./k).*F2)-(((Siga.*P)./k).*F2)-((W./k).*Fa2); %bending moment
dy = ya+((Siga./k).*F2)+((Ma./P).*F3)+((Ra./(k.*P)).*F4)-((W./(k.*P)).*Fa4); %deflection

Q = 2/3*(ro^3-ri^3);
Tau = (V*Q)/(I*dt);
Sig = (M*ro)/I;
Sy = 503e6;

MS = Sy/abs(min(Sig))-1

figure(1)
plot(x,V)
grid on
title('shear force vs. x')
xlabel('x (m)')
ylabel('shear force (N)')
set(gca,'Fontsize', 18)

figure(2)
plot(x,M)
grid on
title('moment vs. x')
xlabel('x (m)')
ylabel('moment (Nm)')
set(gca,'Fontsize', 18)

figure(3)
plot(x,dy)
grid on
title('deflection vs. x')
xlabel('x (m)')
ylabel('deflection (m)')

figure(4)
plot(x,Tau)
grid on
title('Shear stress vs. x')
xlabel('x (m)')
ylabel('Shear stress (Pa)')

figure(5)
plot(x,Sig)
grid on
title('Bending stress vs. x')
xlabel('x (m)')
ylabel('Bending stress (Pa)')