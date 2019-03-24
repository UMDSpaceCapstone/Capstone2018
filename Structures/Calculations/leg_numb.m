clear all

% number of legs
n = 3:10;

%angle of legs
theta = 0:60;

% length of leg (actual length & top view)
l = 1.8;
r = l.*sind(theta);

% Weight of one leg
w1 = 40;

% Weight of n legs
w = w1.*n';

% Angle between legs
ang = 360./n';

% Half angle between legs
ang2 = ang./2;

% Shortest Radius before tipping over
r_crit = (1+r).*cosd(ang2);

% Difference in Shortest Radius 
% before tipping over between designs
for i = 1:length(r_crit(:,1))-1
    dr(i,:) = (r_crit(i,:)./r_crit(i+1,:));
end


%% Plot r_crit vs theta (for 4 legs)

y_crit = 1.7;
y1 = l*cosd(theta(1));

for i = 1:length(theta)
    if y1 > y_crit
        y1 = l*cosd(theta(i));
        theta1 = theta(i);
        r1 = r_crit(2,i);
        j = i;
    end
end

for i = 1:length(theta)
    r11(i) = r1;
end

figure
plot(theta,r_crit(2,:),'.-','MarkerSize',10)
hold on
grid on
plot(theta,r11,'-','MarkerSize',20)
title('Distance to the Edge Vs Leg Angle')
xlabel('Theta (deg)')
ylabel('Edge Distance (m)')
s = ['y = ' num2str(y1)];
legend('Data',s)

%% Plot r_crit vs number of legs

figure
plot(n(2:end)-.5,dr(:,j)','.-','MarkerSize',10)
hold on
grid on
hold off
title('Distance to the Edge Vs Number of Legs')
xlabel('Number of Legs')
ylabel('Edge Distance (m)')

%% 
% 
% d = 1;
% 
% h = (r1-2.9.*sind(theta)+d*cosd(theta))./sind(theta);
% figure
% plot(theta,h)
% title('Height of  Vs Leg Angle')
% xlabel('Theta (deg)')
% ylabel('Edge Distance (m)')