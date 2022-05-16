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
           xs_k = cos(theta(k+1));
           ys_k = sin(theta(k+1));
           xs_k1 = cos(theta(k+2));
           ys_k1 = sin(theta(k+2));
           gamma_k = -2*sin(theta(k+1));
           gamma_k1 = -2*sin(theta(k+2));
           [fa,fb]= panelinf(xs_k,ys_k,xs_k1,ys_k1,xm(i,j),ym(i,j));
           psi(i,j) = psi(i,j)+gamma_k*fa+gamma_k1*fb;
        end      
    end
end

figure
c = -1.75:0.25:1.75;
contour(xm,ym,psi,c);
