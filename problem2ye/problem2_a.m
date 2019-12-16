clc;
clear;
close all;
addpath('functions_p2');
%%



% measuremebt noise
V=diag([1,1]*DEG2RAD).^2;
x_est0=[15;23];
invM=diag([0^2,0^2]);


%%
[x_true,t,z] = truetraj();
N=length(x_true);


figure(1);
f1=plot(x_true(2,1:N),x_true(1,1:N),'r','LineWidth',2);
axis equal
xlabel('y');ylabel('x');
figure(2);
subplot(2,1,1);
plot(z(1,:));
subplot(2,1,2);
plot(z(2,:));

%% Newton-Rhapson

x_est=zeros(2,N);
x_error=zeros(2,N);
cov_list=zeros(2,N);


x_est_temp=x_est0;
invP=invM;
waypoint_total=[];
for i=1:N
    
    [x_est_temp,P,waypoint]=NewtonRhapson(x_est_temp,invP,z(:,i),V)
    invP=inv(P);
    
    x_est(:,i)=x_est_temp;
    x_error(:,i)=x_est_temp-x_true(:,i);
    cov_list(:,i)=diag(P);
    waypoint_total=[waypoint_total,waypoint];
    
    [evector,eval]=eig(P);
    eli=elipse(evector,3.*sqrt(eval),x_est_temp);
    
    figure(1); hold on;
%     f1=scatter(x_est_temp(2),x_est_temp(1),'b');hold on;
%     f2=plot(x_true(2,1),x_true(1,1),'r*');
%     f3=plot(eli(2,:),eli(1,:),'b');
    
end

f4=plot(waypoint_total(2,1),waypoint_total(1,1),'k^');
f5=plot(waypoint_total(2,:),waypoint_total(1,:),'k');
legend([f1,f4,f5],'true','start','est','Location','Best');xlabel('y(km)');ylabel('x(km)');
title('2(a) Newton Rhapson');
axis equal



figure(3);
subplot(2,1,1);
plot(abs(x_error(1,:)));hold on;
plot(3*sqrt(cov_list(1,:)));
xlabel('sec');ylabel('km');
legend('error','3\sigma');
subplot(2,1,2);
plot(abs(x_error(2,:)));hold on;
plot(3*sqrt(cov_list(2,:)));
xlabel('sec');ylabel('km');
sgtitle('2(b) Newton Rhapson');