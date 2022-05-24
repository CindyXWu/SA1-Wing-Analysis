clc
clear all

n = 101;
x = linspace(0,1,n)';
ue = ones(n,1);
gradient = -0.25;
Re = [10^3, 10^4, 10^5];
% First dimension is different x/L values
% Second dimension is different Re values
theta = zeros(n,3);
f = zeros(n,3);
m = zeros(n,3);
% Local Reynolds number based on momentum thickness
Re_theta = zeros(n,3);

% Location of natural transition
int = zeros(3,1);
transition_Re = zeros(3,1);
% Location of laminar separation
ils = zeros(3,1);
separation_Re = zeros(3,1);

% Add constant offset of 1 to ue since ue/U=1 when x/L=0
ue = ue+x*gradient;

% Calculate dimensionless momentum thickness values
% There are 3 such cases, for 3 different (length-based) Re values
% k loops through Re
% i loops through distance along plate x/L
for k=1:3
    for i=2:n
        f(i,k) = f(i-1,k)+ueintbit(x(i-1),ue(i-1),x(i),ue(i));
    end
    theta(:,k) = sqrt((0.45/Re(k))*(ue.^(-6)).*f(:,k));
    m(:,k) = -Re(k)*theta(:,k).^2*gradient; 
    Re_theta(:,k) = Re(k)*ue.*theta(:,k);
end


% Calculate H
H = arrayfun(@thwaites_lookup,m);

% Calculate He
He = arrayfun(@laminar_He,H);

% Find transition and separation point

separation_index = zeros(3,1);
for k=1:3
    laminar = true;
    for i=1:n
        % Transition
        if log(Re_theta(i,k)) >= 18.4*He(i,k)-21.74
            laminar = false;
            disp([x(i) Re_theta(i,k)/1000])
            int(k) = i;
            transition_Re(k) = Re_theta(i,k);
            break
        elseif m(i,k) >= 0.09
            laminar = false;
            ils(k) = i
            separation_Re(k) = Re_theta(i,k);
            break
        end
    end
end

for i=1:3
    if int(i) ~= 0
        disp(['Natural transition at x/L value ' num2str(x(int(i))) ' with Re_theta ' num2str(transition_Re(i))])
    end
    if ils(i) ~= 0
        disp(['Laminar separation at x/L value ' num2str(x(ils(i))) ' with Re_theta ' num2str(separation_Re(i))])
    end 
end