clc;
clear;
close all;
addpath('functions_p2');
%%


invM=diag([0^2,0^2]);
% measuremebt noise
V=diag([1,1]*DEG2RAD).^2;
x_est_temp0=[10;10];


%%
[xTrue,t,z] = truetraj();
N=1000;%length(xTrue);


figure(1);
plot(xTrue(2,1:N),xTrue(1,1:N),'r','LineWidth',2);
axis equal
xlabel('y');ylabel('x');
figure(2);
subplot(2,1,1);
plot(z(1,:));
subplot(2,1,2);
plot(z(2,:));



%% random

invP=diag([0^2,0^2]);
x_est=zeros(2,N);
x_est_temp=estXfromZ(z(:,1));
for i=1:N
    
    
    k=0;
%     while 1
        k=k+1;
        %measurement update
        H=Jacob_h(x_est_temp);
        P=inv(invP+H'*inv(V)*H);
        invP=inv(P);
        K=P*H'*inv(V);
        x_est_temp=x_est_temp+K*(z(:,i)-h(x_est_temp));
        
%         if(  norm(K*(z(:,i)-h(x_est_temp))) )<0.0000000001
%             k
%             break;
%         end
%     end

    x_est(:,i)=x_est_temp;
end

figure(1);hold on;
for i=1:N
    scatter(x_est(2,i),x_est(1,i),'b');hold on;
end
axis equal;