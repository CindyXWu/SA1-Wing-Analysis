clear
clc
close all
Gamma = 1;
nx = 51;
ny = 41;

xrange = linspace(-2.5,2.5,nx)
yrange = linspace(-2,2,ny)
xc = 0.5;
yc = 0.25;
psi = zeros(nx,ny);

for i=1:51
    for j=1:41
       psi(i,j)= psipv(xc,yc,Gamma,xrange(i),yrange(j));
    end
end

c = -0.4:0.2:1.2
[x,y] = meshgrid(xrange,yrange)
figure
mesh(x,y,psi)
hold on
contour(x,y,psi,c)
hold off
grid on
