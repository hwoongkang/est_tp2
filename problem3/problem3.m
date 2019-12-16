clear all;
close all;
clc;
addpath('functions_p1');
problem1_a

%% //Measurement
[true,time,meas] = truetraj();
%% debug
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
%% debug


%% Initial Condition -(a)
A = [0 0; 0 0];         % Continuous System Matrix
Ts = 1;                 % Sampling Time
F = eye(2) + A*Ts;      % Discretized System Matrix
G = [1;1];              % Noise modelling
W = 0.0005;
d2r = pi/180;
V = [1 0;0 1] * d2r * d2r;
omega = wgn(2,2,W);
x_bar = x_est;          % use result of prob1
%% for Plot
x_est_plot = zeros(2,size(time,2)*2);
K_plot = zeros(4,size(time,2));
M_plot = zeros(2,size(time,2)*2);
M = P;
for i = 1:size(time,2)
    %% time update
    x_bar = F * x_bar + G;
    M = F*M*F' + G*W*G';
    
    x_est_plot(:,2*i-1) = x_bar;
    M_plot(:,2*i-1) = diag(M);
        
    %% Measurement update
    K = M*Jacob_h(x_bar)'/V;
    M = inv(inv(M) + Jacob_h(x_est)'/V*Jacob_h(x_est));
    x_bar = x_bar + K*(meas(:,i)-Jacob_h(x_bar)*x_bar);
    
    
    K_plot(:,i) = [K(1) K(2) K(3) K(4)]';
    x_est_plot(:,2*i) = x_bar;
    M_plot(:,2*i) = diag(M);
end


%% debug
time_resize = zeros(1,size(time,2));
k = 0;
for i = 1:size(time,2)*2
   if mod(i,2) == 1
      k=k+1;
      value = time(1,k); 
   end
   time_resize(1,i) = value;
end
figure()
subplot(2,1,1); plot(time_resize,x_est_plot(1,1:end),'k');
title('x_est1'); grid on; xlabel('time(sec)'); ylabel('position(m)'); 
subplot(2,1,2); plot(time_resize,x_est_plot(2,1:end),'k');
title('x_est2'); grid on; xlabel('time(sec)'); ylabel('position(m)');
figure()
titleh = ["K_1","K_2","K_3","K_4"];
for i = 1:4
   subplot(2,2,i); plot(time,K_plot(i,:),'k');
   title(titleh(i)); grid on; xlabel('time(sec)'); ylabel('gain');
end
figure()
titleh = ["P_11","P_22"];
for i =1:2
   subplot(2,1,i); plot(time_resize,M_plot(i,:),'k'); ylabel('RMSE');  xlabel('time(sec)');
end

%% 