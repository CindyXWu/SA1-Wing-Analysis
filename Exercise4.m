clc 
clear all

np = 100;
xmin = 0;
xmax = 5;
ymin = 0;
ymax = 4;
nx=51;
ny=41;
xm = zeros(nx,ny);
ym = zeros(nx,ny);
theta = (0:np)*2*pi/np;
xs = zeros(np);
ys = zeros(np);
gamma = zeros(np);
psi = zeros(nx,ny);

for i=1:nx
    for j=1:ny
        xm(i,j)= xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j)= ymin + (j-1)*(ymax-ymin)/(ny-1);
        for k=0:np-1
           %k is loop over all panels
           xs(k+1) = cos(theta(k+1));
           ys(k+1) = sin(theta(k+1));
           xs(k+2) = cos(theta(k+2));
           ys(k+2) = sin(theta(k+2));
           gamma_k = -2*sin(theta(k+1));
           gamma_k1 = -2*sin(theta(k+2));
           [fa,fb]= panelinf(xs(k+1),ys(k+1),xs(k+2),ys(k+2),xm(i,j),ym(i,j));
           psi(i,j) = psi(i,j)+gamma_k*fa+gamma_k1*fb;
        end      
        %Now add contribution from free stream
        psi(i,j) = psi(i,j)+ym(i,j);
    end
end

figure
hold on
c = -1.75:0.25:1.75;
contour(xm,ym,psi,c);
plot(xs,ys);
hold off

figure
c = -1.75:0.25:1.75;
contour(xm,ym,psi,c);
