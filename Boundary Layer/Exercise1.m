x = linspace(0,1,101);
ue = ones(101,1);
Re=1000;
f = zeros(101,1);
for i=2:101
    f(i) = f(i-1)+(ue(i))^(-6)*ueintbit(x(i-1),ue(i-1),x(i),ue(i));
end

theta = sqrt((0.45/Re)*f);

blasius = 0.664*sqrt(x)/(Re^(0.5));
f_int = sum(f);

figure
hold on
plot(x, theta, 'r--');
plot(x, blasius);
legend('Thwaites','Blasius')
xlabel('x/L')
ylabel('theta/L')
