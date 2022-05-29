function [infa infb] = panelinf(xa,ya,xb,yb,x,y)
    tangential = [xb yb] - [xa ya];
    t = tangential/norm(tangential);
    normal = [(ya-yb) (xb-xa)];
    n = normal/norm(normal);
    r = [x y] - [xa ya];
    X = dot(r,t);
    Y = dot(r,n);
    if abs(Y)<10^(-7);
        Y = 10^(-7);
    end
    [infa infb] = refpaninf(norm(tangential),X,Y);
end


    