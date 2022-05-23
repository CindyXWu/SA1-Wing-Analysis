clc
clear all

n = 101;
% In all following code, we work with the 1st dimension or column dimension
% as n (linspace creates row vectors, so transpose)
x = linspace(0,1,n)';
% 101x3 because each column is for a different gradient
ue = ones(n,1);
gradient = -0.25;
Re = [10^3, 10^4, 10^5];
% First dimension is different x/L values
% Second dimension is different gradients
% Third dimension is different Re values
theta = zeros(n,3);
f = zeros(n,3);
m = zeros(n,3);
% Local Reynolds number based on momentum thickness
Re_theta = zeros(n,3);
% Final 3*3 matrices to output transition values
transition_loc = zeros(3,1);
transition_Re = zeros(3,1);
separation_Re = zeros(3,1);
% Location of turbulent transition
int = zeros(3,1);
% Location of laminar transition
ils = zeros(3,1);

% Fill ue
% Each column of ue is for a different value of gradient wrt dimensionless
% distance

% Add constant offset of 1 to ue since ue/U=1 when x/L=0
ue = ue+x*gradient;

% Calculate dimensionless momentum thickness values
% There are 9 such cases, for 3 different ue values and 3 different
% (length-based) Re values
% k loops through Re
% i loops through distance along plate x/L
for k=1:3
    for i=2:n
        f(i,k) = f(i-1,k)+(ue(i))^(-6)*ueintbit(x(i-1),ue(i-1),x(i),ue(i));
    end
    theta(:,k) = sqrt( (0.45/Re(k)) * f(:,k) );
    m(:,k) = -Re(k)*theta(:,k).^2*gradient;
    
    for i=1:n
        Re_theta(i,k) = Re(k)*ue(i)*theta(i,k);
    end
end


% Calculate H
H = arrayfun(@thwaites_lookup,m);

% Calculate He
He = arrayfun(@laminar_He,H);

% Find transition and separation point
laminar = true;
separation_index = zeros(3,1);
for k=1:3
    for i=1:n
        % Separation
        if m(i,k) >= 0.09
            laminar = false;
            separation_index(k)=i;
            break
        end
    end
    for i=1:n
        % Transition
        if log(Re_theta(i,k)) >= 18.4*He(i,k)-21.74
            laminar = false;
%             disp([x(i) Re_theta(i,k)/1000])
            transition_loc(k) = x(i);
            transition_Re(k) = Re_theta(i,k);
            break
        end
    end
    % Does separation or transition occur first?
    if separation_index(k) < transition_loc(k)
        ils(k) = x(separation_index(k));
        separation_Re = Re(separation_index(k),k);
    else
        int(k) = x(separation_index(k));
    end
end


disp([transition_loc])
disp([transition_Re])
disp([int,ils])
disp([separation_Re])