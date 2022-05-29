function [infa infb] = refpaninf(Delta,X,Y)
    if abs(Y)<10^(-7)
        Y = 10^(-7);
    end
    I0vals = -((X*log(X^2+Y^2)-(X-Delta)*log((X-Delta)^2+Y^2)-2*Delta+2*Y*(atan(X/Y)-atan((X-Delta)/Y))))/(4*pi);
    I1vals = 1/(8*pi)*((X^2+Y^2)*log(X^2+Y^2)-((X-Delta)^2+Y^2)*log((X-Delta)^2+Y^2)-2*X*Delta+Delta^2);
    infa = (1-X/Delta)*I0vals-I1vals/Delta;
    infb = X/Delta*I0vals+I1vals/Delta;
end