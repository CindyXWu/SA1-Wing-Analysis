function rhsvec = build_rhs(xs,ys,alpha)
    %xs is a vector of the edge of panel coordinates
    np = length(xs)-1;
%     psifs = zeros(np+1,np+1);
    %b vector that stores 100X1 values
    rhsvec = zeros(np+1,1);
    for i=1:np+1
        for j=1:np+1
            i1=mod(i,np+1);
            if i1==0
               i1=np+1; 
            end
            j1=mod(i,np+1);
            if j1==0
               j1=np+1; 
            end
            psifs=ys(j)*cos(alpha)-xs(i)*sin(alpha);
            psifs1=ys(j1)*cos(alpha)-xs(i1)*sin(alpha);
            rhsvec(i,1) = psifs-psifs1;
        end
    end
end