close all; clear all;
addpath('functions')

%% Measurement
[true,time,meas] = truetraj();
% debug
figure();
plot(true(2,:),true(1,:),'k');
title('True Trajectory & measurement');
hold on;
plot(meas(2,:),meas(1,:),'b');
grid on;
legend('True','Measurement'); xlabel('y(km)'); ylabel('x(km)'); axis equal;  axis([0 50 0 40]);
% debug


%% Initial Condition
A = [0 0; 0 0];         % Continuous System Matrix
Ts = 1;                 % Sampling Time
F = eye(2) + A*Ts;      % Deiscretized System Matrix
G = [1;1];
% omega = wgn(2,2,W);

