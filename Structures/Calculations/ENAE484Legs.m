% Juliette Abbonizio Lander Legs Calculations
% using the mass and the delta vs calculate the forces
% assuming the legs will be made out of 6061 aluminum

% masses (kg)
mpr_des = .10*22745;
mi_des = 2242.4; %includes structural mass
mtot_asc = 13015;
mtot = mpr_des + mtot_asc + mi_des;

%moon gravity (m/s^2)
g = 1.62;

% velocity at touchdown (m/s)
v = 2;

%impulse (N-s)
I = mtot*v;

%weight
W = mtot*g;

%energy
Ke = 1/2*mtot*v^2;

% choose values for the inner and outer radius (m) (subject to change)
ro = .3
ri = .2

% use area of an annulus (ring)
A = pi*(ro^2 - ri^2);

%moment of intertia of a ring (kg-m^2)
%density = 2700; % kg/m^3
%volume = pi*L*(ro^2-ri^2);
m = 500/4 %(kg)
I = 1/2*m*(ro^2 + ri^2);

%E for 6061 aluminum (Pa)
E = 68.9e9;

% yield strength (sigma) (Pa)
%O = Pcr/A
O = 276e6;
Pcr = O*A

%buckling (using aluminum)
%Pcr = (pi^2*E*I)/(KL)^2
K = .5; % (fixed-fixed)
L = sqrt(((pi^2)*E*I)/(Pcr*K^2))
