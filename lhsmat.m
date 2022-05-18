function lhsmat = build_lhs(xs,ys)
%xs is a vector of the edge of panel coordinates
    np = length(xs)-1;
    psip = zeros(np,np+1);
    %2x2 matrix that stores fa and fb for the (j-1)th and jth panels
    lhsmat = zeros(np+1,np+1);
    
    %build psi matrix
    %i loops over the panel edges at which psi is evaluated. There are np
    %such points

    %j loops over the end points of the panels which produce the
    %streamfunctions
    for i=1:np
        for j = 1:np+1
            if j==1
                [fa fb] = panelinf(xs(j),ys(j),xs(j+1),ys(j+1),xs(i),ys(i));
                psip(i,j)=fa;
            elseif j==np+1
                [faprev fbprev] = panelinf(xs(j-1),ys(j-1),xs(j),ys(j),xs(i),ys(i));
                psip(i,j)=fbprev;
            else
                [fa fb] = panelinf(xs(j),ys(j),xs(j+1),ys(j+1),xs(i),ys(i));
                [faprev fbprev] = panelinf(xs(j-1),ys(j-1),xs(j),ys(j),xs(i),ys(i));
                psip(i,j)=fa+fbprev;
            end
        end
    end
    
    %build A
    for i=1:np+1
        for j=1:np+1
            if i==1
                if j==1
                    lhsmat(i,j)=1;
                else
                    lhsmat(i,j)=0;
                end
            elseif i==np+1
                if j==np+1
                    lhsmat(i,j)=1;
                else
                    lhsmat(i,j)=1;
                end
            else
                lhsmat(i,j)=psi(i+1,:)-psi(i,:);
            
end