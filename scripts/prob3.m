close all; clear all;
addpath('functions')

%% //Measurement
[true,time,meas] = truetraj();
%% //debug
figure();
plot(true(2,:),true(1,:),'k');
title('True Trajectory');
xlabel('y(km)'); ylabel('x(km)'); axis equal;  axis([0 50 0 40]);
figure();
subplot(2,1,1)
plot(meas(1,:),'k');
title('psi1'); xlabel('sample'); ylabel('rad'); 
grid on;
subplot(2,1,2)
plot(meas(2,:),'k');
title('psi2'); xlabel('sample'); ylabel('rad'); 
grid on;
sgtitle('measuremet');
%% //debug


%% //Initial Condition
A = [0 0; 0 0];         % Continuous System Matrix
Ts = 1;                 % Sampling Time
F = eye(2) + A*Ts;      % Discretized System Matrix
G = [1;1];              % Noise modelling
W = 0.005;
omega = wgn(2,2,W);
xbar = [30; 30];

%% // estimation value
xhat = xbar;
/* M =  number 1 value*/

%% //time update
xbar_next = F * xhat + G;


%% //Measurement update


 