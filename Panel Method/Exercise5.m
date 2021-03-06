clear all
close all
clear 

np = 100;
xs = zeros(np+1,1);
ys = zeros(np+1,1);
theta = (0:np)*2*pi/np;
alpha = pi/18;

% Fill xs,ys
for k=1:np+1
    xs(k) = cos(theta(k));
    ys(k) = sin(theta(k));
end

% Build A, b and find gamma using matrix inversion
A = build_lhs(xs,ys);
b = build_rhs(xs,ys,alpha);
gam = A\b;

% Plot of gamma against theta/pi
theta_plot = theta/pi;
figure
plot(theta_plot,gam);
axis([0 2 -2.5 2.5])

%=======================================================
% Produce plot of streamlines as well

xmin = -5;
xmax = 5;
ymin = -4;
ymax = 4;
nx=101;
ny=81;

% xm and ym are the grid points to evaluate psi at
xm = zeros(nx,ny);
ym = zeros(nx,ny);
psi = zeros(nx,ny);

% gammas are evaluated on the surface of the cylinder at panel
% intersections

for i=1:nx
    for j=1:ny
        xm(i,j)= xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j)= ymin + (j-1)*(ymax-ymin)/(ny-1);
        for k=1:np
           %k is loop over all panels
           xs(k) = cos(theta(k));
           ys(k) = sin(theta(k));
           xs(k+1) = cos(theta(k+1));
           ys(k+1) = sin(theta(k+1));
           [fa,fb]= panelinf(xs(k),ys(k),xs(k+1),ys(k+1),xm(i,j),ym(i,j));
           psi(i,j) = psi(i,j)+gam(k)*fa+gam(k+1)*fb;
        end      
        %Now add contribution from free stream
        psi(i,j) = psi(i,j)+ym(i,j)*cos(alpha)-xm(i,j)*sin(alpha);
    end
end

figure
hold on
c = -1.75:0.25:1.75;
contour(xm,ym,psi,c);
plot(xs,ys);
hold off

%=============================================================
% Total circulation = integral of gam over cylinder
% I will just multiply each gam value by the panel length and sum them up

% Panel length
r = sqrt((xs(1)-xs(2))^2+(ys(1)-ys(2))^2);
xgrid = 0:r:np*r;
% Don't need to worry about np+1 index as it's 0 anyway
Gamma = trapz(xgrid,gam);

display(Gamma)


