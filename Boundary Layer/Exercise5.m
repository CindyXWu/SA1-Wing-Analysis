clc
clear all

% Re is overall Reynolds number
% ue0 is free-stream velocity at beginning of integration range
% duedx is the velocity gradient
global Re ue0 duedx
Re = 1e8;
ue0 = 1;
duedx = -0.5;

x0 = 0.01;
n=101;
thick0(1) = 0.023*x0*(Re*x0)^(-1/6);
thick0(2) = 1.83*thick0(1);
[delx thickhist] = ode45(@thickdash, [0 0.99], thick0);

uevec = ones(length(delx),1);
f = zeros(length(delx),1);
% delx is the difference between x/L values and the initial x/L value
x = x0 + delx;

% Power law momentum thickness estimates
theta7 = 0.037*x.*(Re.*x).^(-1/5);
theta9 = 0.023*x.*(Re.*x).^(-1/6);

% Vector of He values
He = thickhist(:,2)./thickhist(:,1);

figure
hold on
plot(x, thickhist(:,1), 'r--','Displayname','PDE solution momentum thickness');
plot(x, thickhist(:,2),'Displayname','PDE solution energy thickness');
legend('Location','northwest')
hold off
title('Exercise 5 thickness plots Re=10^7, gradient=-0.5')
xlabel('x/L')
ylabel('thickness/L')

% figure
% plot(x,He, 'Displayname','He values')
% title('Exercise 5 He plot for separation')
% xlabel('x/L')
% ylabel('thickness/L')