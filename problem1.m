clc;
clear;
close all;
addpath('functions');

N=500;
% x0 ~ N(x_mean,M)
x_mean=[20;30]; %km
M=diag([2^2,2^2]);
% measuremebt noise
V=diag([1,1]*DEG2RAD).^2;


%% simulation
x_true=repmat(x_mean,1,N);
Zmea=zeros(2,N);
for i=1:N
    Zmea(:,i)=h(x_true(:,i))+sqrt(V)*randn(2,1);
end

%% 1(a)
x_est=zeros(2,N);
for i=1:N
    x_est_temp=x_mean;
    k=1;
    while(true)
        
        Z_est=h(x_est_temp);
        J(k)=1/2*(x_est_temp-x_mean)'*inv(M)*(x_est_temp-x_mean)+1/2*(Zmea(:,i)-Z_est)'*inv(V)*(Zmea(:,i)-Z_est);
        k=k+1;
        H_est=Jacob_h(x_est_temp);
        P=inv(inv(M)+H_est'*inv(V)*H_est);
        GRD=inv(M)*(x_est_temp-x_mean)-H_est'*inv(V)*(Zmea(:,i)-Z_est);
        if(norm(GRD)<0.0000000001)
            break;
        end
        x_est_temp=x_est_temp-P*GRD;
    end
    [evector,eval]=eig(P)
    x_est(:,i)=x_est_temp;
end
figure(2);
for i=1:N
    scatter(x_est(1,i),x_est(2,i));hold on;
end
axis equal;

% [evector,eval]=eig(P);
% sqrt(diag(eval))
% 3*sqrt(diag(eval))
% figure(1)
% plot(J);