function psivsheet = psixysheet(GammaA,GammaB,X,Y,Delta)
    I0vals = ((X*log(X^2+Y^2)-(X-Delta)*log((X-Delta)^2+Y^2)-2*Delta+2*Y*(atan(X/Y)-atan((X-Delta)/Y))/(4*pi);
    I1vals = 1/(8*pi)*((X^2+Y^2)*log(X^2+Y^2)-((X-Delta)^2+Y^2)*log((X-Delta)^2+Y^2)-2*X*Delta+Delta^2);
    fa = (1-X/Delta)*I0vals-I1vals/Delta;
    fb = (X/Delta*I0vals+I1vals/Delta;
    psivsheet = GammaA*fa+GammaB*fb;