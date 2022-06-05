function [int ils itr its delstar theta] = bl_solv(x,cp)
    
    % Define global variables
    global Re duedx ue0
    % Re_list=[1e4 1e5 1e6];
    % Redefine the number of panels
    np = length(x);
    
    % Initialise variables at each numbered point in Fig.8
    ue(1) = 0;
    ue(2:np+1) = (1-cp(1:np)).^0.5;
    x_panel(1) = 0;
    x_panel(2:np+1) = x(1:np);
    theta(1:np) = 0;
    delstar(1:np) = 0;
    He = 0;
    f = 0;
    
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

    %% Laminar Flow
    % Find transition and separation point
    laminar = true;
    
    % Location index
    i = 1;
    
    % Start iteration process
    while laminar && i < np+1
        % Update velocity gradient
        duedx = (ue(i+1)-ue(i))/(x_panel(i+1)-x_panel(i));
        f = f + ueintbit(x_panel(i), ue(i), x_panel(i+1), ue(i+1));
        theta(i) =  sqrt((0.45/Re)*(ue(i+1)^(-6))*f);
        Re_theta = Re*ue(i+1)*theta(i);
        % Iterate He
        m = -Re*(theta(i))^2*duedx;
        H = thwaites_lookup(m);
        % Hence update delstar using H and theta
        delstar(i) = H*theta(i);
        He = laminar_He(H);
        
        % After initialising the integral iteration process, we test for
        % natural transition and/or separation
        if log(Re_theta) >= 18.4*He-21.74 % Empirical
            laminar = false;
            int = i;  
        % Test for laminar seperation
        elseif m >= 0.09
            laminar = false;
            ils = i;
            % Set He to laminar separation value
            He = 1.51509;
        end
        i = i + 1;
    end
    %=================================================================
    %% Turbulent Flow for transition or separation
    % Update the delstar array once it leaves laminar regime
    % He below should be that at laminar separation - 1.51509
    delta_pl = He*theta(i-1); % pl means post-laminar, which take the most updated index i results from the previous while loop
    % Under no turbulent separation condition
    while its==0 && i < np+1
        % Calcuate thick vector (theta and delta)
        thick0(1) = theta(i-1);
        thick0(2) = delta_pl;
        ue0 = ue(i-1);
        duedx = (ue(i+1)-ue(i))/(x_panel(i+1)-x_panel(i));
        [delx thickhist] = ode45(@thickdash,[0,x_panel(i+1)-x_panel(i)],thick0);
        
        % For testing of turbulent reattachment, we extract information
        % from several variables from the 'beginning'
        theta(i) = thickhist(end,1);
        delta_pl = thickhist(end,2);
        He = delta_pl/theta(i);
        H = (11*He+15)/(48*He-59); % Empirical
        delstar(i) = H*theta(i);
        % Now, we test for turbulent reattachment (under the condition that laminar separation has occured and reattachment has not taken place)
        if ils ~= 0 && itr == 0 && He > 1.58 % Emprical for He > 1.58
            itr = i;
        % Test further for turbulent separation
        elseif He < 1.46 % Empirical for He < 1.46
            its = i;
            % Given that He and H remains unchanged after turbulent
            % separation
            H = 2.803;
            % Now we can update the full array of theta and delstar
            for k=i:np
                theta(k) = theta(k-1)*(ue(k)/ue(k+1))^(H+2);
                delstar(k) = H*theta(k);
            end
        end
        i = i + 1;
    end
end