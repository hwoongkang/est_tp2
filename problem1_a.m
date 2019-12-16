clc;
clear;
close all;
addpath('functions_p1');

N=500;
% x0 ~ N(x_mean,M)
x_mean=[20;30]; %km
invM=diag([0^2,0^2]);
% measuremebt noise
V=diag([1,1]*DEG2RAD).^2;
x_est_temp0=[10;10];


%% simulation
x_true=repmat(x_mean,1,N);
Zmea=zeros(2,N);
for i=1:N
    Zmea(:,i)=h(x_true(:,i))+sqrt(V)*[cos(2*pi*i/N);sin(2*pi*i/N)];
end

%% 1 sigma


x_est=zeros(2,N);
for i=1:N
    x_est_temp=x_est_temp0;
%     x_est_temp=estXfromZ(Zmea(:,i));
    k=0;
    while(true)
        
        Z_est=h(x_est_temp);
%         J(k)=1/2*(x_est_temp-x_mean)'*invM*(x_est_temp-x_mean)+1/2*(Zmea(:,i)-Z_est)'*inv(V)*(Zmea(:,i)-Z_est);
        k=k+1;
        H_est=Jacob_h(x_est_temp);
        P=inv(invM+H_est'*inv(V)*H_est);
        GRD=invM*(x_est_temp-x_mean)-H_est'*inv(V)*(Zmea(:,i)-Z_est);
        if(norm(GRD)<0.0000000001)
            k
            break;
        end
        x_est_temp=x_est_temp-P*GRD;
    end

    x_est(:,i)=x_est_temp;
end

figure(1);

plot(x_est(1,:),x_est(2,:),'r','LineWidth',2);hold on;
axis equal;

%% simulation
x_true=repmat(x_mean,1,N);
Zmea=zeros(2,N);
for i=1:N
    Zmea(:,i)=h(x_true(:,i))+3*sqrt(V)*[cos(2*pi*i/N);sin(2*pi*i/N)];
end

%% 3 sigma


x_est=zeros(2,N);
for i=1:N
    x_est_temp=x_est_temp0;
%     x_est_temp=estXfromZ(Zmea(:,i));
    k=0;
    while(true)
        
        Z_est=h(x_est_temp);
%         J(k)=1/2*(x_est_temp-x_mean)'*invM*(x_est_temp-x_mean)+1/2*(Zmea(:,i)-Z_est)'*inv(V)*(Zmea(:,i)-Z_est);
        k=k+1;
        H_est=Jacob_h(x_est_temp);
        P=inv(invM+H_est'*inv(V)*H_est);
        GRD=invM*(x_est_temp-x_mean)-H_est'*inv(V)*(Zmea(:,i)-Z_est);
        if(norm(GRD)<0.0000000001)
            k
            break;
        end
        x_est_temp=x_est_temp-P*GRD;
    end

    x_est(:,i)=x_est_temp;
end

figure(1);

plot(x_est(1,:),x_est(2,:),'g','LineWidth',2);hold on;
axis equal;




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
    while(true)
        
        Z_est=h(x_est_temp);
%         J(k)=1/2*(x_est_temp-x_mean)'*invM*(x_est_temp-x_mean)+1/2*(Zmea(:,i)-Z_est)'*inv(V)*(Zmea(:,i)-Z_est);
        k=k+1;
        H_est=Jacob_h(x_est_temp);
        P=inv(invM+H_est'*inv(V)*H_est);
        GRD=invM*(x_est_temp-x_mean)-H_est'*inv(V)*(Zmea(:,i)-Z_est);
        if(norm(GRD)<0.0000000001)
            k
            break;
        end
        x_est_temp=x_est_temp-P*GRD;
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