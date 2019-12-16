clc;
clear;
close all;
addpath('functions_p1');

N=3;

% measuremebt noise
V=diag([1,1]*DEG2RAD).^2;
% x0 ~ (x_est0,inf);
x_est0=[15;23];
invM=diag([0^2,0^2]);
%% simulation
x_true=repmat([20;30],1,N);
Zmea=zeros(2,N);
for i=1:N
    Zmea(:,i)=h(x_true(:,i))+sqrt(V)*randn(2,1);
end
%% Newton-Rhapson

x_est=zeros(2,N);
for i=1:N
    
    [x_est,P,waypoint]=NewtonRhapson(x_est0,invM,Zmea(:,i),V)
    
    
    
    [evector,eval]=eig(P);
    eli=elipse(evector,3.*eval,x_est);
    
    
    
    figure(1); hold on;
    scatter(x_est(2),x_est(1),'b');hold on;
    plot(x_true(2,1),x_true(1,1),'r*');
    legend('est','true');
    plot(eli(2,:),eli(1,:),'b');
    plot(waypoint(2,:),waypoint(1,:),'k');
    axis equal
    legend('est','true','3\sigma','Newton-Rhapson','Location','Best');
    xlabel('y(km)');ylabel('x(km)');
end



