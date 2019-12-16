clc;
clear;
close all;
addpath('functions_p1');

N=50;
% measuremebt noise
V=diag([1,1]*DEG2RAD).^2;
% x0 ~ (x_est0,inf);
x_est_before_NR=[15;23];
invM_before_NR=diag([0^2,0^2]);
%% simulation
x_true=repmat([20;30],1,N);

Zmea=zeros(2,N);
for i=1:N
    Zmea(:,i)=h(x_true(:,i))+sqrt(V)*randn(2,1);
end
%% x_est0 ±¸ÇÏ±â

[x_est_after_NR,M_after_NR]=NewtonRhapson(x_est_before_NR,invM_before_NR,Zmea(:,1),V)



%% Kalman filter


x_est=zeros(2,N-1);
x_error=zeros(2,N-1);
cov_list=zeros(2,N-1);


x_est_temp=x_est_after_NR;
P=M_after_NR;

waypoint=x_est_temp;
for i=1:N-1
    
    
    %measurement update
    H=Jacob_h(x_est_temp);
    P=inv(inv(P)+H'*inv(V)*H);
    K=P*H'*inv(V);
    x_est_temp=x_est_temp+K*(Zmea(:,i+1)-h(x_est_temp));    
    
    x_est(:,i)=x_est_temp;
    x_error(:,i)=x_est(:,i)-x_true(:,i);
    cov_list(:,i)=diag(P);
    waypoint=[waypoint,x_est_temp];
    
    [evector,eval]=eig(P);
    eli=elipse(evector,3.*sqrt(eval),x_est_temp);
    
    
    figure(1); hold on;
    f1=scatter(x_est_temp(2),x_est_temp(1),'b');hold on;
    f2=plot(x_true(2,1),x_true(1,1),'r*');
    f3=plot(eli(2,:),eli(1,:),'b');
    
end

f4=plot(waypoint(2,1),waypoint(1,1),'k^');
f5=plot(waypoint(2,:),waypoint(1,:),'k');
legend([f1,f2,f3,f4,f5],'est','true','3\sigma','start','Kalman Filter','Location','Best');xlabel('y(km)');ylabel('x(km)');
title('1(b) Kalman Filter');
axis equal



figure(2);
subplot(2,1,1);
plot(abs(x_error(1,:)));hold on;
plot(3*sqrt(cov_list(1,:)));
xlabel('sec');ylabel('km');
legend('error','3\sigma');
subplot(2,1,2);
plot(abs(x_error(2,:)));hold on;
plot(3*sqrt(cov_list(2,:)));
xlabel('sec');ylabel('km');
sgtitle('1(b) Kalman Filter');


