function [int ils itr its delstar theta] = bl_solv(x,cp)
    
    % Define global variables
    global Re duedx ue0
    % Re_list=[1e4 1e5 1e6];
    % Redefine the number of panels
    np = length(x);
    
    % Initialise variables at each numbered point in Fig.8
    ue(1) = 0;
    ue(2:np+1) = (1-cp(1:np)).^0.5;
    x_panel(1:np+1) = 0;
    x_panel(2:np+1) = x(1:np);
    
%     delstar(1:np) = 0;
%     He = 0;
    duedx = zeros(np,1);
    
    % Initialise integral, theta and m
    % theta(1:np) = 0;
    f = zeros(np,1);
    % m = zeros(np,1);
    
    % Initialise indices tracker
    % Location of laminar transition
    int = 0;
    transition_Re = 0;
    % Location of laminar separation
    ils = 0;
    separation_Re = 0;
    % Location of turbulent reattachment
    itr = 0;
    turbulent_reattachment_Re = 0;
    % Location of turbulent separation
    its = 0;
    turbulent_separation_Re = 0;


    % Calculate dimensionless momentum thickness values
    % i loops through distance along plate x/L
    for i=2:np
        duedx(i) = (ue(i)-ue(i-1))/(x_panel(i)-x_panel(i-1));
        f(i) = f(i-1) + ueintbit(x_panel(i-1),ue(i-1),x_panel(i),ue(i));
    end
    theta = sqrt((0.45/Re).*(ue.^(-6)).*f);
    m = -Re*(theta.^2).*duedx;
    Re_theta = Re*ue.*theta;
    % Calculate H
    H = arrayfun(@thwaites_lookup,m);   
    % Calculate delstar
    delstar = H.*theta;
    % Calculate He
    He = arrayfun(@laminar_He,H);
    % Set arbitrary value for first point of He (step i)
    He(1) = 1.57258;
    %=================================================
    %% Laminar Flow
    % Find transition and separation point
    laminar = true;
    for i=1:np
        % Transition
        if log(Re_theta(i)) >= 18.4*He(i)-21.74
            laminar = false;
            int = i;
%             transition_Re = Re_theta(i);
            break
        % Separation
        elseif m(i) >= 0.09
            % If laminar separation detected, He fixed at its separation
            % value (step ii)
            He(i:end)=1.51509;
%             delta = He.*theta;
            laminar = false;
            ils = i;
%             separation_Re = Re_theta(i);
            break  
        end
    end

    %=================================================================
    %% Turbulent Flow
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
    while(its == 0 & i < 101)
        [delx thickhist] = ode45(@thickdash, [0,x(i)-x(i-1)], thick0);
        thick0 = thickhist(end,:);
        duedx = (ue(i+1)-ue(i))/(x_panel(i+1) - x_panel(i));
        % m_i = -Re*(thick0(1))^2*duedx;
        % Fill in theta array for plotting later
        theta(i) = thick0(1);
        % Fill in He array for plotting later
        He(i) = thick0(2)/thick0(1);
        H(i) = (11*He(i) + 15)/(48*He(i) -59);
        % calculate theta and hence delstar
        theta(i) = thickhist(end,1);
        delstar(i) = H(i)*theta(i);

        % Test for turbulent re-attachment
        % Condition: laminar separation has already occurred, and He > 1.58
        % If turbulent re-attachment has already occurred, then do not test
        % again
        if (ils ~= 0 & itr == 0 & He_i >= 1.58)
            itr = i;
%             turbulent_reattachment_Re = Re_theta(i);
        end

        % Test for turbulent separation regardless of if turbulent
        % re-attachment has occurred
        if (He_i <= 1.46)
            its = i;
%             turbulent_separation_Re = Re_theta(i);
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
        % H remains unchanged after turbulent separation
        H = 2.803;
        while k <= 100
            theta(k+1) = theta(k)*(ue(k)/ue(k+1))^(H+2);
            delstar(k) = H*theta(k);
            k = k+1;
        end
    end
end