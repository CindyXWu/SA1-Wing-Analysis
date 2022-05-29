clc
clear global
clear all
close

global Re duedx ue0
Re_list=[1e4 1e5 1e6];
stylelist = ['g','k','r--'];
for j=1:3
    n = 101;
    x = linspace(0,1,n)';
    ue = ones(n,1);
    
    duedx = -0.25;
    Re = Re_list(j);
    % First dimension is different x/L values
    theta = zeros(n,1);
    f = zeros(n,1);
    m = zeros(n,1);
    % Local Reynolds number based on momentum thickness
    Re_theta = zeros(n,1);
    delta = zeros(n,1);

    % Location of natural transition
    int = 0;
    transition_Re = 0;
    % Location of laminar separation
    ils = 0;
    separation_Re = 0;

    % Initialise turbulent trackers (step i)
    % Location of turbulent reattachment
    itr = 0;
    turbulent_reattachment_Re = 0;
    % Location of turbulent separation
    its = 0;
    turbulent_separation_Re = 0;

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

    % Set arbitrary value for first point of He (step i)
    He(1) = 1.57258;

    %=================================================

    % Find transition and separation point
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
            % If laminar separation detected, He fixed at its separation
            % value (step ii)
            He(i:end)=1.51509;
            delta = He.*theta;
            laminar = false;
            ils = i;
            separation_Re = Re_theta(i);
            break  
        end
    end
    
    %=================================================================

    % Turbulent section!
    % Set index at which the ODE solver following laminar separation should follow
    % on. Depends on if separation or transition occurs first.
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

    % Provide initial values for turbulent calculations (step ii)
    theta_sep = theta(i-1);
    % delta = He*theta (step iii)
    delta_sep = 1.51509*theta_sep;
    thick0 = [theta_sep, delta_sep];
    ue0 = ue(i-1);
    
    % Calculations for non-laminar regions: while loop proceeds as long as
    % its = 0 and position index < max value
    while(its == 0 & i <= 101)
        [delx thickhist] = ode45(@thickdash, [0,x(i)-x(i-1)], thick0);
        thick0 = thickhist(end,:);
        He_i = thick0(2)/thick0(1);
        m_i = -Re*(thick0(1))^2*duedx;
        % Fill in theta array for plotting later
        theta(i) = thick0(1);
        % Fill in He array for plotting later
        He(i) = He_i;

        % Test for turbulent re-attachment
        % Condition: laminar separation has already occurred, and He > 1.58
        % If turbulent re-attachment has already occurred, then do not test
        % again
        if (ils ~= 0 & itr == 0 & He_i >= 1.58)
            itr = i;
            turbulent_reattachment_Re = Re_theta(i);
        end

        % Test for turbulent separation regardless of if turbulent
        % re-attachment has occurred
        if (He_i <= 1.46)
            its = i;
            turbulent_separation_Re = Re_theta(i);
            % H unchanged from separation value
            H(its:end) = thwaites_lookup(m(its));
            % Populate He with constant values (separated)
            He(its:end) = He_i
        end
        ue0 = ue(i);
        i = i+1;
    end
    
    % Calculate theta after turbulent separation
    if its ~= 0
        k = its;
        % Calculation for after turbulent separation
        while k <= 100
            theta_a = theta(k);
            uea =  ue(k);
            ueb = ue(k+1);
            theta_b = theta_a*(uea/ueb)^(H(k)+2);
            theta(k+1) = theta_b;
            k = k+1;
        end
    end

%     hold on
%     plot(x, theta,stylelist(j),'Displayname',strcat('Re=',num2str(Re)));
%     legend('Location','northwest')
%     xlabel('x/L')
%     title('Exercise 6 momentum thickness plot gradient=0')
%     ylabel('theta/L')
    
    hold on
    plot(x, He, stylelist(j),'Displayname',strcat('Re=',num2str(Re)));
    xlabel('x/L')
    title('Exercise 6 He plot gradient=0')
    ylabel('He')
    legend('Location','northwest')
    
    if int ~= 0
        disp(['Natural transition at x/L value ' num2str(x(int)) ' with Re_theta ' num2str(transition_Re) ' and ReL=' num2str(Re)])
    end
    if ils ~= 0
        disp(['Laminar separation at x/L value ' num2str(x(ils)) ' with Re_theta ' num2str(separation_Re) ' and ReL=' num2str(Re)])
    end
    if itr ~= 0
        disp(['Turbulent reattachment at x/L value ' num2str(x(itr)) ' with Re_theta' num2str(turbulent_separation_Re) ' and ReL=' num2str(Re)])
    end
    if its ~= 0
        disp(['Turbulent separation at x/L value ' num2str(x(its)) ' with Re_theta' num2str(turbulent_reattachment_Re) ' and ReL=' num2str(Re)]) 
    end
    disp('-------------------------------------------------------------------------------------------------------------------------------------')
end

% xline(0.49, 'c','Displayname','Natural transition for Re=10^6')
% xline(0.5, 'Color','#A2142F','Displayname','Laminar separation for Re=10^4')
% xline(0.59, 'Color','#EDB120','Displayname','Turbulent re-attachment for Re=10^5')
% xline(0.82, 'm','Displayname','Turbulent separation for Re=10^4')

% xline(0.37, 'm','Displayname','Natural transition for Re=10^7')
% For reference: turbulent separation at x/L=1 for Re=10^5 from duedx =
% -0.382