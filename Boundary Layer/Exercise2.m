x = linspace(0,1,101);
ue = ones(101,3);
gradients = [-0.1,0,0.1];
Re = [5*10^6, 10^7, 2*10^7];
% First dimension is different x/L values
% Second dimension is different gradients
% Third dimension is different Re values
f = zeros(101,3,3);
theta = zeros(101,3,3);
m = zeros(101,3,3);
Re_theta = zeros(101,3,3);

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
        for i=2:101
            f(i,j,k) = f(i-1,j,k)+(ue(i,j))^(-6)*ueintbit(x(i-1),ue(i-1,j),x(i),ue(i,j));
        end
        theta(:,j,k) = sqrt((0.45/Re(k))*f);
        m(:,j,k) = -Re(k)*theta(:,j,k).^2*gradients(j);
        
        for i=1:101
             Re_theta = Re(k)*ue(i,j)*
    end
end


% Calculate H
H = arrayfun(thwaites_lookup,m);

% Calculate He
He = arrayfun(laminar_He,m);

blasius = 0.664*sqrt(x)/(Re^(0.5));
f_int = sum(f);

figure
hold on
plot(x, theta, 'r--');
plot(x, blasius);
legend('Thwaites','Blasius')
xlabel('x/L')
ylabel('theta/L')