function [int ils itr its delstar theta] = bl_solv(x,cp)
    
    % Define global variables
    global Re duedx ue0
 
    % Redefine the number of panels
    np = length(x);
    
    % Initialise variables at each numbered point in Fig.8
    ue(1:np) = (1-cp).^0.5;
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
    duedx = (ue(1)-0)/(x(1)-0);

    f = f + ueintbit(0, 0, x(1), ue(1));
    theta(1) =  sqrt((0.45/Re)*(ue(1)^(-6))*f);
    Re_theta = Re*ue(1)*theta(1);
    % Iterate He
    m = -Re*(theta(1))^2*duedx;
    H = thwaites_lookup(m);
    % Hence update delstar using H and theta
    delstar(1) = H*theta(1);
    He = laminar_He(H);
    % Start iteration process
    while laminar && i < np
        i = i + 1;
        % Update velocity gradient
        duedx = (ue(i)-ue(i-0))/(x(i)-x(i-1));
        f = f + ueintbit(x(i-1), ue(i-1), x(i), ue(i));
        theta(i) =  sqrt((0.45/Re)*(ue(i)^(-6))*f);
        Re_theta = Re*ue(i)*theta(i);
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
    end
    %=================================================================
    delta_pl = He*theta(i); % pl means post-laminar, which take the most updated index i results from the previous while loop
  
    %% Turbulent Flow for transition or separation
    % Update the delstar array once it leaves laminar regime
    % He below should be that at laminar separation - 1.51509
   
    % Under no turbulent separation condition
    while its==0 && i < np
        i = i + 1;
        ue0 = ue(i-1);
        % Calcuate thick vector (theta and delta)
        thick0(1) = theta(i-1);
        thick0(2) = delta_pl;
        duedx = (ue(i)-ue(i-1))/(x(i)-x(i-1));
        [delx thickhist] = ode45(@thickdash,[0,x(i)-x(i-1)],thick0);
        thick0 = thickhist(end,:);
        % For testing of turbulent reattachment, we extract information
        % from several variables from the 'beginning'
        theta(i) = thickhist(end,1);
        delta_pl = thickhist(end,2);
        He = delta_pl/theta(i);
        if He >= 1.46
            H = (11*He+15)/(48*He-59);
        else
            H = 2.803;
        end 
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
            for k=i+1:np
                theta(k) = theta(k-1)*(ue(k-1)/ue(k))^(H+2);
                delstar(k) = H*theta(k);
            end
        end
    end
end