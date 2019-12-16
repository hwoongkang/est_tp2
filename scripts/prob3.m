addpath('../functions_p1/functions_p1')

%% Measurement


%% Initial Condition;
A = [0 0; 0 0];         % Continuous System Matrix
Ts = 1;                 % Sampling Time
F = eye(2) + A*Ts;      % Deiscretized System Matrix
G = [1;1];
omega = wgn(2,2,W);
