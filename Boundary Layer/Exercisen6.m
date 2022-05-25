clc
clear all
global Re duedx ue0
Re_list=[1e6 1e7];
stylelist = ['g','k'];
for j=1:2
    n = 101;
    x = linspace(0,1,n)';
    ue = ones(n,1);
    
    duedx = 0;
    Re = Re_list(j);
    % First dimension is different x/L values
    theta = zeros(n,1);
    f = zeros(n,1);
    m = zeros(n,1);
    % Local Reynolds number based on momentum thickness
    Re_theta = zeros(n,1);
    delta = zeros(n,1);
    itr = 0;
    its = 0;

    % Location of natural transition
    int = zeros(3,1);
    transition_Re = zeros(3,1);
    % Location of laminar separation
    ils = zeros(3,1);
    separation_Re = zeros(3,1);

    % Add constant offset of 1 to ue since ue/U=1 when x/L=0
    ue = ue+x*duedx;

    % Calculate dimensionless momentum thickness values
    % i loops through distance along plate x/L
    for i=2:n
        f(i) = f(i-1)+ueintbit(x(i-1),ue(i-1),x(i),ue(i));
    end
    theta = sqrt((0.45/Re)*(ue.^(-6)).*f);
    m = -Re*theta.^2*duedx; 
    Re_theta = Re*ue.*theta;

    % Calculate H
    H = arrayfun(@thwaites_lookup,m);

    % Calculate He
    He = arrayfun(@laminar_He,H);
    He(1) = 1.57258;
    % Find transition and separation point

    separation_index = zeros(3,1);
        laminar = true;
    for i=1:n
        % Transition
        if log(Re_theta(i)) >= 18.4*He(i)-21.74
            laminar = false;
            int = i;
            transition_Re = Re_theta(i);
            break
        % Separation
        elseif m(i) >= 0.09
            He=1.51509;
            delta = He.*theta;
            laminar = false;
            ils = i;
            separation_Re = Re_theta(i);
            break  
        end
    end

    if(int == 0 & ils ~= 0)
    % Laminar separation
        i = ils+1;
    elseif(int ~= 0 & ils == 0)
    % Laminar-turbulent transition
        i = int+1;
    elseif(int ~= 0 & ils ~= 0)
    % Check which occurs first
        i = min([int ils])+1;
    else
        % Fully laminar - set i such that while loop does not run
        i = 102;
    end
    theta_sep = theta(i-1);
    delta_sep = 1.51509*theta_sep;
    thick0 = [theta_sep, delta_sep];
    ue0 = ue(i-1);
    while(its == 0 & i <= 101)
        [delx thickhist] = ode45(@thickdash, [0,x(i)-x(i-1)], thick0);
        thick0 = thickhist(end,:);
        He_i = thick0(2)/thick0(1);
        m_i = -Re*thickhist(1)*duedx;
        % Fill in theta array for plotting later
        theta(i) = thick0(1);
        % Fill in He array for plotting later
        He(i) = He_i;
        % Fill in m array for reference
        m(i) = m_i;
        % Test for laminar separation followed by re-attachment
        if(ils ~= 0 & He_i >= 1.58)
            itr = i;
        end
        % Test for turbulent separation
        if(int ~= 0 & m >= 0.09)
            its = i;
        end
        i = i+1;
        ue0 = ue(i-1);
    end

    if its ~= 0
        i = its;
        % Calculation for after turbulent separation
        while i <= 101
            theta_a = theta(i);
            uea =  ue(i);
            ueb = ue(i+1);
            H = thwaites_lookup(m(i));
            theta_b = theta_a*(uea/ueb)^(H+2);
            theta(i+1) = theta_b;
            i = i+1;
        end
    end

    hold on
    plot(x, theta,stylelist(j),'Displayname',strcat('Re=',num2str(Re)));
    legend('Location','northwest')
    xlabel('x/L')
    title('Exercise 6 momentum thickness plot gradient=0')
    ylabel('theta/L')
% 
%     hold on
%     plot(x, He, stylelist(j),'Displayname',strcat('Re=',num2str(Re)));
%     legend('Location','northwest')
%     xlabel('x/L')
%     title('Exercise 6 He plot gradient=0')
%     ylabel('theta/L')

    if int ~= 0
        disp(['Natural transition at x/L value ' num2str(x(int)) ' with Re_theta ' num2str(transition_Re)])
    end
    if ils ~= 0
        disp(['Laminar separation at x/L value ' num2str(x(ils)) ' with Re_theta ' num2str(separation_Re)])
    end
end