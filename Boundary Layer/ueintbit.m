function f = ueintbit(xa,ua,xb,ub)
    ubar = 0.5*(ua+ub);
    delta_u = ub - ua;
    delta_x = xb - xa;
    f = (ubar^5+5/6*ubar^3*delta_u^2+1/16*ubar*delta_u^4)*delta_x;
end