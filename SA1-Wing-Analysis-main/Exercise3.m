clear
clc
close all
Gamma = 3;
xa=3.5;
ya=2.5;
xb=1.6;
yb=1.1;
xmin = 0;
xmax = 5;
ymin = 0;
ymax = 4;
nv = 100;
nx = 51;
ny = 41;
psi = zeros(nx,ny);
xm = zeros(nx,ny);
ym = zeros(nx,ny);
Delta = 1.5;
fa = zeros(nx,ny);
fb = zeros(nx,ny);

for i=1:nx
    for j=1:ny
        xm(i,j)= xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j)= ymin + (j-1)*(ymax-ymin)/(ny-1);
        [fa(i,j),fb(i,j)]= panelinf(xa,ya,xb,yb,xm(i,j),ym(i,j));
    end
end

c = -0.15:0.05:0.15
figure
contour(xm,ym,fa,c)

figure
contour(xm,ym,fb,c)

%start of discretised
dl = Delta/nv;
gamma = zeros(nv,1);
discrete_fa = zeros(nx,ny); %psi at one location due to many point vortices
discrete_fb = zeros(nx,ny);

for i=1:nx
    for j=1:ny
        xm(i,j)= xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j)= ymin + (j-1)*(ymax-ymin)/(ny-1);
        for k=1:nv
            xc = xa+(xb-xa)*(2*k-1)/(2*nv);
            yc = ya+(yb-ya)*(2*k-1)/(2*nv);
            r = sqrt((yb-ya)^2+(xb-xa)^2);
            gamma_a = r*(nv+0.5-k)/(nv^2);
            gamma_b = r*(k-0.5)/(nv^2);
            discrete_fa(i,j) = discrete_fa(i,j) + psipv(xc,yc,gamma_a,xm(i,j),ym(i,j));
            discrete_fb(i,j) = discrete_fb(i,j) + psipv(xc,yc,gamma_b,xm(i,j),ym(i,j));
        end
    end
end

figure
c = -0.15:0.05:1.15;
contour(xm,ym,discrete_fa,c);

figure
contour(xm,ym,discrete_fb,c);
