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
            f(i,j,k) = f(i-1,j,k)+ ueintbit(x(i-1),ue(i-1,j),x(i),ue(i,j));
        end
        theta(:,j,k) = sqrt( (0.45/Re(k))*((ue(:,j)).^(-6)).*f(:,j,k) );
        m(:,j,k) = -Re(k)*theta(:,j,k).^2*gradients(j);
        
        for i=1:n
            Re_theta(i,j,k) = Re(k)*ue(i,j)*theta(i,j,k);
        end
    end
end

% x = linspace(0,1,n)';
% duedxtest = -0.1;
% ftest = zeros(n,1);
% ue0 = 1;
% uetest = ue0+duedxtest*x;
% Retest = 5e6;
% thetatest = zeros(n,1);
% mtest = zeros(n,1);
% Re_thetatest = zeros(n,1);
% 
% 
% for i=2:n
%     ftest(i) = ftest(i-1)+ueintbit(x(i-1),uetest(i-1),x(i),uetest(i));
% end
% thetatest = sqrt( (0.45/Retest)*(((uetest).^(-6)).*ftest ));
% mtest = -Retest*thetatest.^2*duedxtest;
% 
% for i=1:n
%     Re_thetatest(i) = Retest*uetest(i)*thetatest(i);
% end

% Calculate H
H = arrayfun(@thwaites_lookup,m);

% Calculate He
He = arrayfun(@laminar_He,H);

% Find transition point
for j=1:3
    for k=1:3
        laminar = true;
        for i=1:n
            if log(Re_theta(i,j,k)) >= 18.4*He(i,j,k)-21.74
                disp(i)
                laminar = false;
                disp([x(i) Re_theta(i,j,k)/1000])
                transition_loc(j,k) = x(i);
                transition_Re(j,k) = Re_theta(i,j,k);
                break
            end
        end
    end
end

% laminar = true;
% for i=1:n
%     if log(Re_thetatest(i)) >= 18.4*He(i)-21.74
%         disp(i)
%         laminar = false;
%         disp([x(i) Re_thetatest(i)/1000])
%         transition_loc = x(i);
%         transition_Re = Re_thetatest(i);
%         break
%     end
% end

disp(transition_loc)
disp(transition_Re)