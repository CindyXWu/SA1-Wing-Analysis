function dthickdx = thickdash(xmx0,thick)
    global Re ue0 duedx
    ue = ue0+duedx*xmx0;
    He = thick(1)/thick(2);
    if He >= 1.46
        H = (11*He+15)/(48*He-59);
    else
        H = 2.803;
    end   
    Re_theta = Re*ue*thick(1);
    cf = 0.091416*((H-1)*Re_theta)^(-0.232)*exp(-1.26*H);
    c_diss = 0.010012*((H-1)*Re_theta)^(-1/6);
    dthickdx(1,1) = cf/2-((H+2)./ue).*duedx*thick(1);
    dthickdx(2,1) = c_diss-(3./ue).*duedx*thick(2);
end