function rhsvec = build_rhs(xs,ys,alpha)
    % xs is a vector of the edge of panel coordinates
    % This has np+1 entries as the first and last entries are repeated
    % There are np unique coordinates
    np = length(xs)-1;
    % rhsvec: vector that stores 101X1 values
    rhsvec = zeros(np+1,1);
    % Equations 7 and 8 require 0 on RHS
    rhsvec(1)=0;
    rhsvec(np+1)=0;
    % Equation 6 fills in the rest of the vector
    for i=2:np
        rhsvec(i) = ys(i)*cos(alpha)-xs(i)*sin(alpha)
        -ys(i+1)*cos(alpha)+xs(i+1)*sin(alpha);
    end
end


