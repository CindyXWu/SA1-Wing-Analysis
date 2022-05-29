clc
clear all
clear figure

% Re is overall Reynolds number
% ue0 is free-stream velocity at beginning of integration range
% duedx is the velocity gradient
global Re ue0 duedx
Re = 1e7;
ue0 = 1;
duedx = 0;
x0 = 0.01;
n=101;
thick0(1) = 0.023*x0*(Re*x0)^(-1/6);
thick0(2) = 1.83*thick0(1);
[delx thickhist] = ode45(@thickdash, [0 0.99], thick0);

uevec = ones(length(delx),1);
f = zeros(length(delx),1);
x = x0 + delx;

theta7 = 0.037*x.*(Re.*x).^(-1/5);
theta9 = 0.023*x.*(Re.*x).^(-1/6);

figure
hold on
plot(x, thickhist(:,1), 'r--','Displayname','PDE solution');
plot(x, theta7, 'b+','Displayname','Turbulent 1/7th power law');
plot(x, theta9,'g.','Displayname','Turbulent 1/9th power law');
legend('Location','northwest')
hold off
xlabel('x/L')
title('Exercise 4 momentum thickness plot')
ylabel('theta/L')
