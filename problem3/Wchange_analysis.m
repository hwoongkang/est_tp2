clear all; close all;
addpath('functions_p1');
x_est = [19.805459672505162;  29.912102830854717]; % 1번의 결과
P = [0.262585111631561  -0.086710045763324; -0.086710045763324   0.392070074203004]; % 1번의 결과
%% Measurement
[true,time,meas] = truetraj();

%% Initial Condition -(a)
A = [0 0; 0 0];         % Continuous System Matrix
Ts = 1;                 % Sampling Time
F = eye(2) + A*Ts;      % Discretized System Matrix


% G = [1;1];              % Noise modelling with no disturbability
G = [1;1.001];          % disturbability

w = [0.0001 0.001 0.01 0.1 1];
d2r = pi/180;
V = [1 0;0 1] * d2r * d2r;
monte = size(w,2);

error_plot = cell(monte,1);
x_est_plot = zeros(2,size(time,2));
K_plot = zeros(2,2,size(time,2));
K_Plot = cell(monte,1);
for j = 1:size(w,2)
    M = P;
    W = w(j)^2;
    % van loan 
    C = [-F G*W*G'; zeros(2,2) F'];
    vanloan = expm(C);
    W = vanloan([3 4],[3 4])'*vanloan([1 2],[3 4]);
    x_bar = x_est + [10 0]';          % use result of prob1
    for i = 1:size(time,2)
        %% time update
        x_bar = F * x_bar;
        M = F*M*F' + W;
        
        %% Measurement update
        K = M*Jacob_h(x_bar)'/V;
        M = inv(inv(M) + Jacob_h(x_bar)'/V*Jacob_h(x_bar));
        x_bar = x_bar + K*(meas(:,i)-h(x_bar));
        x_est_plot(:,i) = x_bar;
        K_plot(:,:,i) = K;
    end
    errpos = x_est_plot-true;
    error_plot{j,1} = errpos;
    K_Plot{j,1} = K_plot;
end

%% plot
% kalman gain
title_K = ["K_{11}","K_{12}","K_{21}","K_{22}"];
legend_K = [sprintf("%s%s","W= ",num2str(w(1)^2,'%.2e')),sprintf('%s%s',"W= ",num2str(w(2)^2,'%.2e')),sprintf('%s%s',"W= ",num2str(w(3)^2,'%.2e')),sprintf('%s%s',"W= ",num2str(w(4)^2,'%.2e')),sprintf('%s%s',"W= ",num2str(w(5)^2,'%.2e'))];
figure();
for i =1:monte
    K_11 = reshape(K_Plot{i}(1,1,:),1,size(time,2));
    K_12 = reshape(K_Plot{i}(1,2,:),1,size(time,2));
    K_21 = reshape(K_Plot{i}(2,1,:),1,size(time,2));
    K_22 = reshape(K_Plot{i}(2,2,:),1,size(time,2));
    plot_K = {K_11 K_12 K_21 K_22};
    for j = 1:4
        subplot(2,2,j); plot(time,plot_K{j}); title(title_K(j)); xlabel('time(s)'); ylabel('gain'); grid on; hold on;
        if i == monte
            for k = 1:4
                legend(legend_K(1),legend_K(2),legend_K(3),legend_K(4),legend_K(5));
            end
        end
    end
end
sgtitle('Kalman gain');

% pos error
legend_pos = [sprintf("%s%s","W= ",num2str(w(1)^2,'%.2e')),sprintf('%s%s',"W= ",num2str(w(2)^2,'%.2e')),sprintf('%s%s',"W= ",num2str(w(3)^2,'%.2e')),sprintf('%s%s',"W= ",num2str(w(4)^2,'%.2e')),sprintf('%s%s',"W= ",num2str(w(5)^2,'%.2e'))];
figure();
for i = 1:monte 
    for j =1:2
        poserr = reshape(error_plot{i}(j,:),1,size(time,2));
        subplot(2,1,j); plot(time,poserr); hold on; grid on; xlabel('time(s)'); ylabel('error(km)');
        if i == monte
           legend(legend_pos(1),legend_pos(2),legend_pos(3),legend_pos(4),legend_pos(5));
        end
    end
end
sgtitle('Position error');