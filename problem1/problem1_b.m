clc;
clear;
close all;
addpath('functions_p1');

N=3;
% x0 ~ N(x_mean,M)
x_mean=[20;30]; %km
invM=diag([0^2,0^2]);
% measuremebt noise
V=diag([1,1]*DEG2RAD).^2;

%% simulation
x_true=repmat(x_mean,1,N);
Zmea=zeros(2,N);
for i=1:N
    Zmea(:,i)=h(x_true(:,i))+sqrt(V)*randn(2,1);
end

%% random


x_est=zeros(2,N);
for i=1:N
    x_est_temp=[10;10];
    invP=diag([0^2,0^2]);
    k=0;
    while 1
        k=k+1;
        %measurement update
        H=Jacob_h(x_est_temp);
        P=inv(invP+H'*inv(V)*H);
        K=P*H'*inv(V);
        x_est_temp=x_est_temp+K*(Zmea(:,i)-h(x_est_temp));
        
        if(  norm(K*(Zmea(:,i)-h(x_est_temp))) )<0.0000000001
            k
            break;
        end
    end

    x_est(:,i)=x_est_temp;
end

figure(1);
for i=1:N
    scatter(x_est(1,i),x_est(2,i),'b');hold on;
end
axis equal;

% [evector,eval]=eig(P);
% sqrt(diag(eval))
% 3*sqrt(diag(eval))
% figure(1)
% plot(J);