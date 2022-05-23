clc
clear all

n = 101;
% In all following code, we work with the 1st dimension or column dimension
% as n (linspace creates row vectors, so transpose)
x = linspace(0,1,n)';
% 101x3 because each column is for a different gradient
ue = ones(n,3);
gradients = [-0.1,0,0.1];
Re = [5*10^6, 10^7, 2*10^7];
% First dimension is different x/L values
% Second dimension is different gradients
% Third dimension is different Re values
theta = zeros(n,3,3);
f = zeros(n,3,3);
m = zeros(n,3,3);
% Local Reynolds number based on momentum thickness
Re_theta = zeros(n,3,3);
% Final 3*3 matrices to output transition values
transition_loc = zeros(3);
transition_Re = zeros(3);

% Fill ue
% Each column of ue is for a different value of gradient wrt dimensionless
% distance

for i = 1:3
    % Add constant offset of 1 to ue since ue/U=1 when x/L=0
    ue(:,i) = ue(:,i)+x*gradients(i);
end

% Calculate dimensionless momentum thickness values
% There are 9 such cases, for 3 different ue values and 3 different
% (length-based) Re values
% j loops through gradients
% k loops through Re
% i loops through distance along plate x/L
for j=1:3
    for k=1:3
        for i=2:n
            f(i,j,k) = f(i-1,j,k)+ue(i,j)^(-6)*ueintbit(x(i-1),ue(i-1,j),x(i),ue(i,j));
        end
        theta(:,j,k) = sqrt( (0.45/Re(k)) * f(:,j,k) );
        m(:,j,k) = -Re(k)*theta(:,j,k).^2*gradients(j);
        
        for i=1:n
            Re_theta(i,j,k) = Re(k)*ue(i,j)*theta(i,j,k);
        end
    end
end

% Calculate H
H = arrayfun(@thwaites_lookup,m);

% Calculate He
He = arrayfun(@laminar_He,H);

% Find transition point
laminar = true;
for j=1:3
    for k=1:3
        for i=1:n
            if log(Re_theta(i,j,k)) >= 18.4*He(i,j,k)-21.74
                laminar = false;
                disp([x(i) Re_theta(i,j,k)/1000])
                transition_loc(j,k) = x(i);
                transition_Re(j,k) = Re_theta(i,j,k);
                break
            end
        end
    end
end

disp(transition_loc)
disp(transition_Re)