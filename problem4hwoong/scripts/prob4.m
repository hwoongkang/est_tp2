clc;
clear;
close all;
%% prob 4.
addpath("functions")
addpath("../problem1/functions_p1")

%% discretization
% continuous state matrix
A = [0,1,0,0;
    0,0,0,0;
    0,0,0,1;
    0,0,0,0];

% continuous Gamma Matrix (with Q_c = W) Q_c=W=(km/s^2/sqrt(Hz))^2
gamma = [0;1;0;1];

% Gamma * W * Gamma^T
Q = gamma * gamma.';

% Van Loan's Discretization
C = expm([-A,Q;zeros(4), A.']);

Ad = C(5:8,5:8).';
Qd = Ad * C(1:4,5:8);
%% part (a)
[xTrue,t,z]= truetraj();

R = diag([DEG2RAD, DEG2RAD].^2);

[x0pos,P0pos]=NewtonRhapson([10;10],zeros(2),z(:,1),R);

% allocation
x = zeros(4,length(xTrue));
x([1,3],1) = x0pos;
x([2,4],1) = [50/3600;0];

x0 = x(:,1);

P0 = zeros(4);

P0([1,3],[1,3]) = P0pos;
P0([2,4],[2,4]) = 1E4 * eye(2);

P = zeros(4,4,length(xTrue));
P(:,:,1) = P0;

M = zeros(4,4,length(xTrue));

getXWithW = @(W) getX(Ad,Qd,W,R,xTrue,z,x0,P0);
getXKWithW = @(W) getXK(Ad,Qd,W,R,xTrue,z,x0,P0);
getXWithW1AndW2 = @(W1,W2) getX2(Ad,Qd,W1,W2,R,xTrue,z,x0,P0);
getXKWithW1AndW2 = @(W1,W2) getXK2(Ad,Qd,W1,W2,R,xTrue,z,x0,P0);


%%
err = 1E20;
WAns = 0;
f = waitbar(0,'W...');



for W = [0.001^2,0.01^2, 0.02^2,0.03^2,0.05^2, 0.1^2 ,1^2 ] %[1E-5:1E-5:1E-3] %9E-4%
    waitbar(W,f,sprintf("current best: %.3f^2",WAns));
    [x,K_list] = getXKWithW(W^2);
    poserr = x([1,3],:) - xTrue;
    errnow = trace(poserr.' * poserr) / length(x);
    if errnow<err
        err = errnow;
        WAns = W;
        bestXWithW = x;
    end
    
    figure;sgtitle(sprintf("4(a) W=%.3f^2",sqrt(W)));
    for i=1:4
        for j=1:2
            subplot(4,2,(i-1)*2+j);plot(t(2:end),K_list(:,i,j));title(sprintf("K(%d,%d)",i,j));
        end
    end
    
    figure;hold on; title(sprintf("4(a) Trajectory W=%.3f^2",sqrt(W)));
    plot(x(3,:),x(1,:)); plot(xTrue(2,:),xTrue(1,:));
    legend('est','true');
end
fprintf("4(a) best Answer %f \n",WAns);


close(f)



%%


x = getXWithW1AndW2(0.00025^2,0.0002^2);

g = waitbar(0,'W1&W2...');

W1Ans = 0;
W2Ans = 0;
err2 = 1E20;

% best(W1, W2): (5.18, 2.67) * 1E-4
for W1 = [0.00001^2, 0.0228^2, 0.05^2, 0.1^2   ]%5.17*1E-4:1E-6:5.19*1E-4
    
    waitbar(W1*1000/2,g,sprintf("current best: %f,%f",W1Ans,W2Ans));
    
    for W2 =[0.001^2, 0.0163^2, 0.05^2, 0.1^2]  % 2.6*1E-4:1E-5:2.8*1E-4
        
        [x2,K_list] = getXKWithW1AndW2(W1^2,W2^2);
        poserr2 = x2([1,3],:) - xTrue;
        errnow = trace(poserr2.' * poserr2)/length(x2);
        if errnow<err2
            bextXWithW1AndW2 = x2;
            err2 = errnow;
            W1Ans = W1;
            W2Ans = W2;
        end
        
        figure;sgtitle(sprintf("4(b) Kalman Gain W1=%.5f^2 W2=%.5f^2",sqrt(W1),sqrt(W1)));
        for i=1:4
            for j=1:2
                subplot(4,2,(i-1)*2+j);plot(t(2:end),K_list(:,i,j));title(sprintf("K(%d,%d)",i,j));
            end
        end
        
        figure;hold on; title(sprintf("4(b) Trajectory W1=%.5f^2 W2=%.5f^2",sqrt(W1),sqrt(W2)));
        plot(x(3,:),x(1,:)); plot(xTrue(2,:),xTrue(1,:));
        legend('est','true');
    end
end

fprintf("4(b) best Answer W1=%f W2=%f \n",W1Ans,W2Ans);

close(g)
function x = getX(Ad,Qd,W,R,xTrue,z,x0,P0)
    % tuning parameter
    Qd = W*Qd;
    x = zeros(4,length(xTrue));
    x(:,1) = x0;
    P = zeros(4,4,length(xTrue));
    P(:,:,1) = P0;
    for ind = 2:length(xTrue)
        % time update
        x(:,ind) = Ad * x(:,ind-1);
        M = Ad * P(:,:,ind-1) * Ad.' + Qd;
        
        % measurement update
        dz = z(:,ind) - h(x([1,3],ind));
        Htemp = Jacob_h(x([1,3],ind));
        H = zeros(2,4);
        H(:,[1,3]) = Htemp;
        
        K = M * H.' * inv(H*M* H .' +R);
        P(:,:,ind) = (eye(4) - K*H) * M;
        
        dx = K*dz;
        
        x(:,ind) = x(:,ind) + dx;
    end
end
function [x,K_list] = getXK(Ad,Qd,W,R,xTrue,z,x0,P0)
    % tuning parameter
    Qd = W*Qd;
    x = zeros(4,length(xTrue));
    K_list= zeros(length(xTrue)-1,4,2);
    x(:,1) = x0;
    P = zeros(4,4,length(xTrue));
    P(:,:,1) = P0;
    for ind = 2:length(xTrue)
        % time update
        x(:,ind) = Ad * x(:,ind-1);
        M = Ad * P(:,:,ind-1) * Ad.' + Qd;
        
        % measurement update
        dz = z(:,ind) - h(x([1,3],ind));
        Htemp = Jacob_h(x([1,3],ind));
        H = zeros(2,4);
        H(:,[1,3]) = Htemp;
        
        K = M * H.' * inv(H*M* H .' +R);
        P(:,:,ind) = (eye(4) - K*H) * M;
        
        dx = K*dz;
        
        K_list(ind-1,:,:)=K;
        x(:,ind) = x(:,ind) + dx;
    end
end
function x = getX2(Ad,Qd,W1,W2,R,xTrue,z,x0,P0)
    % van loan again
    A = [0,1,0,0;
    0,0,0,0;
    0,0,0,1;
    0,0,0,0];
    gamma = [0,0;1,0;0,0;0,1];
    W = [W1,0;0,W2];
    C = expm([-A,gamma * W * gamma.';zeros(4),A.']);
    Ad = C(5:8,5:8).';
    Qd = Ad * C(1:4,5:8);
    
    x = zeros(4,length(xTrue));
    x(:,1) = x0;
    P = zeros(4,4,length(xTrue));
    P(:,:,1) = P0;
    for ind = 2:length(xTrue)
        % time update
        x(:,ind) = Ad * x(:,ind-1);
        M = Ad * P(:,:,ind-1) * Ad.' + Qd;
        
        % measurement update
        dz = z(:,ind) - h(x([1,3],ind));
        Htemp = Jacob_h(x([1,3],ind));
        H = zeros(2,4);
        H(:,[1,3]) = Htemp;
        
        K = M * H.' * inv(H*M * H .' +R);
        P(:,:,ind) = (eye(4) - K*H) * M;
        
        dx = K*dz;
        
        x(:,ind) = x(:,ind) + dx;
    end
end
function [x,K_list] = getXK2(Ad,Qd,W1,W2,R,xTrue,z,x0,P0)
    % van loan again
    A = [0,1,0,0;
    0,0,0,0;
    0,0,0,1;
    0,0,0,0];
    gamma = [0,0;1,0;0,0;0,1];
    W = [W1,0;0,W2];
    C = expm([-A,gamma * W * gamma.';zeros(4),A.']);
    Ad = C(5:8,5:8).';
    Qd = Ad * C(1:4,5:8);
    
    x = zeros(4,length(xTrue));
    K_list= zeros(length(xTrue)-1,4,2);
    x(:,1) = x0;
    P = zeros(4,4,length(xTrue));
    P(:,:,1) = P0;
    for ind = 2:length(xTrue)
        % time update
        x(:,ind) = Ad * x(:,ind-1);
        M = Ad * P(:,:,ind-1) * Ad.' + Qd;
        
        % measurement update
        dz = z(:,ind) - h(x([1,3],ind));
        Htemp = Jacob_h(x([1,3],ind));
        H = zeros(2,4);
        H(:,[1,3]) = Htemp;
        
        K = M * H.' * inv(H*M * H .' +R);
        P(:,:,ind) = (eye(4) - K*H) * M;
        
        dx = K*dz;
        
        K_list(ind-1,:,:)=K;
        x(:,ind) = x(:,ind) + dx;
    end
end
