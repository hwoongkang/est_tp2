%% function to provide the true traj and measurements

function [x,t,z] = truetraj(seed)
    if nargin<1
        seed = 1;
    end
    rng(seed);
    
    PHR2PSEC = 1/3600; % km/hr (per hour) to km/s (per sec)
    %% constants
    % Ts = 1
    dt = 1;

    % radius of the trajectory
    R0 = 10;
    
    % origin of the ''
    Origin = [20;30];
    
    % initial velocity
    V0 = PHR2PSEC * 50;    
    
    % initial yaw
    THETA0 = 0;
    
    % initial position
    x0 = Origin + R0 * [cos(THETA0); sin(THETA0)];
    
    % measurement noise std
    nu = 1 * DEG2RAD;
    
    %% variables
    % with some margin
    t = 0:dt:75*60;
        
    % allocation
    x = zeros(2,length(t));
    z = zeros(size(x));
    % variable to hold the vel
    v = V0;
    
    % to hold yaw
    theta = THETA0;
    
    x(:,1) = x0;
    z(:,1) = h(x0) + nu *randn(2,1);
    %% true traj 
    % iteration!
    for ind = 2:1:length(t)
        % get acc profile
        a = acc(t(ind));
        
        % euler integration
        v = v + a* dt;
        
        % convert it to angular velocity,
        % to get a perfect circle trajectory
        omega = v / R0;
        
        theta = theta - omega * dt;
        
        % and update the pos
        x(:,ind) = Origin + R0 * [cos(theta);sin(theta)];
        % measurement
        z(:,ind) = h(x(:,ind)) + nu*randn(2,1);
        % if returned
        if (ind>64*60)&& (theta<-2*pi)
            break
        end
    end
    % crop
    t = t(1:ind);
    x = x(:,1:ind);
    z = z(:,1:ind);
end