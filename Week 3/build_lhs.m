function lhsmat = build_lhs(xs,ys)
% xs is a vector of the edge of panel coordinates
    np = length(xs)-1;
    psip = zeros(np,np+1);
    lhsmat = zeros(np+1,np+1);
    
    % Build psi matrix
    % i loops over the panel edges at which psi is evaluated. There are np
    % such points

    % j loops over the end points of the panels which produce the
    % streamfunctions
    for i = 1:np  
        fb = 0;
        for j = 1:np
            psip(i,j) = fb;
            [fa, fb] = panelinf(xs(j),ys(j),xs(j+1),ys(j+1),xs(i),ys(i));
            psip(i,j)=fa+psip(i,j);
        end
        psip(i,np+1) = fb;
    end
   
    
    % Build A
    % A has a dimension of 101x101, i.e. np+1xnp+1
%     lhsmat(1,1)=1;
    % multiply psip by gamma for zero incidence (theoretical values) to see
    % if it gives constant psi values on RHS
    % Excluding top and bottom rows, for later modification
    for i=1:np-1
        lhsmat(i+1,:)=psip(i+1,:)-psip(i,:);
    end
%     lhsmat(np+1,np+1)=1;
    
    % First row
    lhsmat(1,1) = 1;
    lhsmat(1,2) = -1;
    lhsmat(1,3) = 0.5;
    lhsmat(1,np-1)= -0.5;
    lhsmat(1,np)= 1;

    % np+1 row
    lhsmat(np+1,2) = 1;
    lhsmat(np+1,3) = -0.5;
    lhsmat(np+1,np-1)= 0.5;
    lhsmat(np+1,np)= -1;
    lhsmat(np+1,np+1) = 1;
end