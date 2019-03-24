% time (s)
dt = [0:.0001:5];

% Mass (kg)
m = 1.77476e4;
% Damping constant (kg/s)
c = 4*0.35*2*m;
% Spring constant (N/m)
k = 710.764e3;
% Gravitational constant (m/s^2)
g = 1.62;
% length of spring (m)
l = 0.4;

y(1) = 0;
v(1) = 1;


for i = 2:length(dt)
 a(i-1) = (-c/m)*v(i-1)-(k/m)*y(i-1)+4*g;
 y(i) = y(i-1)+v(i-1)*.0001;
 v(i) = v(i-1)+a(i-1)*.0001;
end

y = l-y;

plot(dt,y(1:length(dt)))

xlabel('t(s)')
ylabel('y(m)')