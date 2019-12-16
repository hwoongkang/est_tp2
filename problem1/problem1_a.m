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
x_est_temp0=[15;23];


%% simulation
x_true=repmat(x_mean,1,N);
Zmea=zeros(2,N);
for i=1:N
    Zmea(:,i)=h(x_true(:,i))+sqrt(V)*randn(2,1);
end
%% random

x_est=zeros(2,N);
for i=1:N
    x_est_temp=x_est_temp0;
%     x_est_temp=estXfromZ(Zmea(:,i));
    k=0;
    waypoint=x_est_temp;
    while(true)
        
        Z_est=h(x_est_temp);
%         J(k)=1/2*(x_est_temp-x_mean)'*invM*(x_est_temp-x_mean)+1/2*(Zmea(:,i)-Z_est)'*inv(V)*(Zmea(:,i)-Z_est);
        k=k+1;
        H_est=Jacob_h(x_est_temp);
        P=inv(invM+H_est'*inv(V)*H_est);
        GRD=invM*(x_est_temp-x_mean)-H_est'*inv(V)*(Zmea(:,i)-Z_est);
        norm(GRD)
        if(norm(GRD)<0.0000000001)
            k
            break;
        end
        x_est_temp=x_est_temp-P*GRD;
        waypoint=[waypoint,x_est_temp]; 
    end

    x_est(:,i)=x_est_temp;
    [evector,eval]=eig(P);

    
     eli=elipse(evector,3.*eval,x_est(:,i));
    figure(1); hold on;
    scatter(x_est(2,i),x_est(1,i),'b');hold on;
    plot(x_true(2,:),x_true(1,:),'r*');
    legend('est','true');
    plot(eli(2,:),eli(1,:),'b');
    plot(waypoint(2,:),waypoint(1,:),'k');
    axis equal
    legend('est','true','3\sigma','Newton-Rhapson','Location','Best');
    xlabel('y(km)');ylabel('x(km)');
end

