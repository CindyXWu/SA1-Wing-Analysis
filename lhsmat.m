function lhsmat = build_lhs(xs,ys)
%xs is a vector of the edge of panel coordinates
    np = length(xs)-1;
    psip = zeros(np,np+1);
    %2x2 matrix that stores fa and fb for the (j-1)th and jth panels
    lhsmat = zeros(np+1,np+1);
    %Checked
    
    
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
    %Checked
    
    %build A
    %A has a dimension of 100x100, ie np+1xnp+1
    for i=1:np+1
        for j=1:np+1
            %prevent index exceeding range -> when i=np+1+1, it should go
            %back to 1
            i1=mod(i,np+1);
            if i1==0
               i1=np+1; 
            end
            lhsmat(i1,j)=psip(i1+1,j)-psip(i1,j);
        end
    end
end