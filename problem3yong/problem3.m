clear all;
close all;
clc;
addpath('functions_p1');
problem1_a

%% Measurement
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
%% 


%% Initial Condition -(a)
A = [0 0; 0 0];         % Continuous System Matrix
Ts = 1;                 % Sampling Time
F = eye(2) + A*Ts;      % Discretized System Matrix


G = [1;1];              % Noise modelling
w = 0.002:0.0001:0.01;
d2r = pi/180;
V = [1 0;0 1] * d2r * d2r;


%% compare monte
f = waitbar(0,'processing...');
error_plot = zeros(1,size(w,2));
x_est_plot = zeros(2,size(time,2));
for j = 1:size(w,2)
    M = P;
    W = w(j)^2;
    x_bar = x_est + [10 0]';          % use result of prob1
    for i = 1:size(time,2)
        %% time update
        x_bar = F * x_bar;
        M = F*M*F' + G*W*G';
        
        %% Measurement update
        K = M*Jacob_h(x_bar)'/V;
        M = inv(inv(M) + Jacob_h(x_bar)'/V*Jacob_h(x_bar));
        x_bar = x_bar + K*(meas(:,i)-h(x_bar));
        x_est_plot(:,i) = x_bar;
    end
    waitbar(j/size(w,2),f,'prob3-(a) finding w...');
    rootmeansq = rms(x_est_plot-true,2);
    error_plot(:,j) = sqrt(rootmeansq(1)^2 + rootmeansq(2)^2);
end
close(f);
[value, index] = min(error_plot);
W_1 = w(index)^2;                                                 % rmse가  최소가 되는 w 선정
figure();
plot(w.^2,error_plot); ylabel('RMSE(km)'); xlabel('W1');
titterton = sprintf("%s%s",'Minimum RMSE at W = ',num2str(W_1,'%.6f'));
title(titterton);
%% for Plot
W=W_1;
x_est_plot = zeros(2,size(time,2)*2);
K_plot = zeros(4,size(time,2));
M_plot = zeros(2,size(time,2)*2);
meas_residual = zeros(1,size(time,2));
M = P;
x_bar = x_est + [10 0]';
for i = 1:size(time,2)
    %% time update
    x_bar = F * x_bar;
    M = F*M*F' + G*W*G';
    
    x_est_plot(:,2*i-1) = x_bar;
    M_plot(:,2*i-1) = diag(M);
        
    %% Measurement update
    K = M*Jacob_h(x_bar)'/V;
    M = inv(inv(M) + Jacob_h(x_bar)'/V*Jacob_h(x_bar));
    meas_residual(:,i) = norm(meas(:,i)-h(x_bar));
    x_bar = x_bar + K*(meas(:,i)-h(x_bar));
    K_plot(:,i) = [K(1) K(2) K(3) K(4)]';
    x_est_plot(:,2*i) = x_bar;
    M_plot(:,2*i) = diag(M);
end


%% plot
time_resize = zeros(1,size(time,2)*2);
true_resize = zeros(2,size(time,2)*2);
k = 0;
for i = 1:size(time,2)*2
   if mod(i,2) == 1
      k=k+1;
      value = time(1,k); 
      value2 = true(:,k);
   end
   time_resize(1,i) = value;
   true_resize(:,i) = value2;
end
figure()
plot(x_est_plot(2,:),x_est_plot(1,:),'b','Linewidth',0.5);
hold on;
plot(true(2,:),true(1,:),'k');
title('Trajectory'); grid on; xlabel('time(sec)'); ylabel('position(m)');
legendary = sprintf('%s%s','W = ',num2str(W,'%.6f'));
legend(legendary);
figure()
titleh = ["K_1","K_2","K_3","K_4"];
for i = 1:4
   subplot(2,2,i); plot(time,K_plot(i,:),'r');
   title(titleh(i)); grid on; xlabel('time(sec)'); ylabel('gain');
   axis([0 4500 -30 30]);
end
sgtitle('Kalman gain');
figure()
titleh = ["P_{11}","P_{22}"];
for i =1:2
   subplot(2,1,i); plot(time_resize,x_est_plot(i,:)-true_resize(i,:),'r'); hold on; plot(time_resize,3*sqrt(M_plot(i,:)),'k'); hold on; plot(time_resize,3*-sqrt(M_plot(i,:)),'k'); ylabel('distance(km)');  xlabel('time(sec)'); 
   legend('error','covariance(3\sigma)');
   title(titleh(i));
end
sgtitle('Covariance');
figure()
plot(time,meas_residual);

title('Measurement residual(optional)');
ylabel('rad'); xlabel('time(sec)');

%% Initial Condition -(b)
A = [0 0; 0 0];         % Continuous System Matrix
Ts = 1;                 % Sampling Time
F = eye(2) + A*Ts;      % Discretized System Matrix
G = [1 0;0 1];
w1 = 0.04:0.005:0.1;     % trial-error 로 이쁜구간 찾음
w2 = 0.04:0.005:0.1;     % trial-error 로 이쁜구간 찾음
d2r = pi/180;
V = [1 0;0 1] * d2r * d2r;



%% compare monte
f = waitbar(0,'processing');
error_plot2 = zeros(1,size(W,2));
x_est_plot = zeros(2,size(time,2));
Z = ones(size(w1,2),size(w2,2));
for k = 1:size(w1,2)
    W1 = w1(k)^2;
    for j = 1:size(w2,2)
        M = P;
        W2 = w2(j)^2;
        W = [W1 0;0 W2];
        x_bar = x_est + [10 0]';          % use result of prob1
        for i = 1:size(time,2)
            %% time update
            x_bar = F * x_bar;
            M = F*M*F' + G*W*G';

            %% Measurement update
            K = M*Jacob_h(x_bar)'/V;
            M = inv(inv(M) + Jacob_h(x_bar)'/V*Jacob_h(x_bar));
            x_bar = x_bar + K*(meas(:,i)-h(x_bar));
            x_est_plot(:,i) = x_bar;
        end
        rootmeansq = rms(x_est_plot-true,2);
        error_plot2(:,size(w1,2)*(k-1)+j) = sqrt(rootmeansq(1)^2 + rootmeansq(2)^2);
        waitbar((size(w1,2)*(k-1)+j)/(size(w1,2)*size(w2,2)),f,'prob 3-(b). finding w1 and w2...');
    end
end
close(f)
figure()
[X,Y] = meshgrid(w1.^2,w2.^2);
for ind_tmp = 1:size(w1,2)
    for ind_tmp2 = 1:size(w2,2)
       Z(ind_tmp,ind_tmp2) = error_plot2(size(w1,2)*(ind_tmp-1)+ind_tmp2);
    end
end

surf(X,Y,Z);
xlabel('W1');
ylabel('W2');
zlabel('RMSE(km)');
[value, index] = min(error_plot2);
k = floor(index/size(w1,2)) + 1;
if k==size(w1,2)+1
   k = size(w1,w); 
end
j = mod(index,size(w1,2));
if j==0
   j = 8; 
end
W1 = w1(k)^2;
W2 = w2(j)^2;
hwang = sprintf('%s%s%s%s','Minimum RMSE at W1= ',num2str(W1),' and W2=',num2str(W2));
title(hwang);
W = [W1 0; 0 W2];
%% for Plot
W1=0.0098278^2;
W2=0.0098278^2;
W = [W1 0; 0 W2];
x_est_plot = zeros(2,size(time,2)*2);
K_plot = zeros(4,size(time,2));
M_plot = zeros(2,size(time,2)*2);
meas_residual = zeros(1,size(time,2));
M = P;
x_bar = x_est + [10 0]';
for i = 1:size(time,2)
    %% time update
    x_bar = F * x_bar;
    M = F*M*F' + G*W*G';
    
    x_est_plot(:,2*i-1) = x_bar;
    M_plot(:,2*i-1) = diag(M);
        
    %% Measurement update
    K = M*Jacob_h(x_bar)'/V;
    M = inv(inv(M) + Jacob_h(x_bar)'/V*Jacob_h(x_bar));
    meas_residual(:,i) = norm(meas(:,i)-h(x_bar));
    x_bar = x_bar + K*(meas(:,i)-h(x_bar));
    K_plot(:,i) = [K(1) K(2) K(3) K(4)]';
    x_est_plot(:,2*i) = x_bar;
    M_plot(:,2*i) = diag(M);
end

%% Plot
figure()
plot(x_est_plot(2,:),x_est_plot(1,:),'b','Linewidth',0.5);
hold on;
plot(true(2,:),true(1,:),'k');
title('trajectory'); grid on; xlabel('time(sec)'); ylabel('position(m)');
legendary1 = sprintf('%s%s','W1 = ',num2str(W1,'%.6f'));
legendary2 = sprintf('%s%s','W2 = ',num2str(W2,'%.6f'));
legend(sprintf('%s%s%s',legendary1,", ",legendary2));
figure()
titleh = ["K_1","K_2","K_3","K_4"];
for i = 1:4
   subplot(2,2,i); plot(time,K_plot(i,:),'r');
   title(titleh(i)); grid on; xlabel('time(sec)'); ylabel('gain');
   axis([0 4500 -30 30]);
end
sgtitle('Kalman gain');
figure()
titleh = ["P_11","P_22"];
for i =1:2
   subplot(2,1,i);  plot(time_resize,x_est_plot(i,:)-true_resize(i,:),'r'); hold on;  plot(time_resize,3*sqrt(M_plot(i,:)),'k'); hold on; plot(time_resize,3*-sqrt(M_plot(i,:)),'k'); ylabel('distance(km)');  xlabel('time(sec)'); 
   legend('error','covariance(3\sigma)');
   title(titleh(i));
end
sgtitle('Covariance');
figure()
plot(time,meas_residual);
title('Measurement residual(optional)');
ylabel('rad'); xlabel('time(sec)');