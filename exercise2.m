clear
clc
close all
Gamma = 3;
xmin = -2.5;
xmax = 2.5;
ymin = -2;
ymax = 2;
nv = 100;
nx = nv;
ny = nv;
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
        [fa(i,j),fb(i,j)]= refpaninf(Delta,xm(i,j),ym(i,j));
    end
end

c = -0.15:0.05:0.15;
figure
contour(xm,ym,fa,c)

figure
contour(xm,ym,fb,c)

dl = Delta/nv;
gamma = zeros(nv,1);
yc = 0;
discrete_fa = zeros(nv); %psi at one location due to many point vortices
discrete_fb = zeros(nv);


for i=1:nx
   
    for j=1:ny
        xm(i,j)= xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j)= ymin + (j-1)*(ymax-ymin)/(ny-1);
        for k=1:nv
            xc = dl*(2*k-1)/2;
            gamma_a = dl*(nv+0.5-k)/nv;
            gamma_b = dl*(k-0.5)/nv;
   
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

