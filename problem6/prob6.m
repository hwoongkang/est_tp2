clear all; close all; clc;
addpath('functions_p1');
addpath('functions');
%% discretization
% continuous state matrix
A = [0,1,0,0;
    0,0,0,0;
    0,0,0,1;
    0,0,0,0];

% continuous Gamma Matrix (with Q_c = W)
gamma = [0;1;0;1];

% Gamma * W * Gamma^T
Q = gamma * gamma.';

% Van Loan's Discretization
C = expm([-A,Q;zeros(4), A.']);

Ad = C(5:8,5:8).';
Qd = Ad * C(1:4,5:8);

R = diag([DEG2RAD, DEG2RAD].^2);
[xTrue,t,z]= truetraj();

[x0pos,P0pos]=NewtonRhapson([10;10],zeros(2),z(:,1),R);

%% state 10 น่, covariance 10น่
x0pos_10 = xTrue(:,1) + 10*(x0pos - xTrue(:,1));
x([1,3],1) = x0pos_10;
x([2,4],1) = [0;0];
x0 = x(:,1);
P0pos_10 = P0pos * 10^2;
P0 = zeros(4);
P0([1,3],[1,3]) = P0pos_10;
P0([2,4],[2,4]) = 1E4 * eye(2);
P = zeros(4,4,length(xTrue));
P(:,:,1) = P0;


%% (b)
getXWithW1AndW2 = @(W1,W2) getX2(A,Qd,W1,W2,R,xTrue,z,x0,P0);

M = zeros(4,4,length(xTrue));
% initial error var
err = 1E20;
WAns = 0;

g = waitbar(0,'W1&W2...');

W1Ans = 0;
W2Ans = 0;
err2 = 1E20;

% best(W1, W2): (5.175, 2.67) * 1E-4
for W1 = 5.175*1E-4 %5.17*1E-4:1E-7:5.19*1E-4
    
    waitbar(W1*1000/2,g,sprintf("current best: %f,%f",W1Ans,W2Ans));
    
    for W2 = 2.67*1E-4 %2.6*1E-4:1E-6:2.8*1E-4
        
        % mean square error
        [x2,K_plot] = getXWithW1AndW2(W1^2,W2^2);
        poserr2 = x2([1,3],:) - xTrue;
        errnow = trace(poserr2.' * poserr2)/length(x2);
        if errnow<err2
            bextXWithW1AndW2 = x2;
            err2 = errnow;
            W1Ans = W1;
            W2Ans = W2;
        end
    end
end

close(g)

bestX2 = bextXWithW1AndW2;

figure
plot(bestX2(3,:),bestX2(1,:),'.k')
hold on
plot(xTrue(2,:),xTrue(1,:),'--r')
title(sprintf("(b) best result with W1 = %.3e, W2 = %.3e",W1Ans,W2Ans))
axis equal

figure()
subplot(2,1,1); plot(t,poserr2(1,:),'k'); title('Error of position_x'); xlabel('time(s)'); ylabel('error(km)'); grid on;
subplot(2,1,2); plot(t,poserr2(2,:),'k'); title('Error of position_y'); xlabel('time(s)'); ylabel('error(km)'); grid on;
sgtitle('position error');

figure()
titleh = ["K_{11}","K_{12}","K_{21}","K_{22}","K_{31}","K_{32}","K_{41}","K_{42}"];
for i = 1:8
    tmpmat = reshape(K_plot(1,1,:),1,4065);
    subplot(4,2,i); plot(t(50:end),tmpmat(50:end)); title(titleh(i)); xalbel('time(s)'); ylabel('gain'); grid on;
end
sgtitle('Kalman gain(from 50seconds)');

%% functions
function x = getX(Ad,Qd,W,R,xTrue,z,x0,P0)
    %Qd
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

function [x,K_plot] = getX2(A,Qd,W1,W2,R,xTrue,z,x0,P0)
    % tuning parameter
    % van loan again
    gamma = [0,0;1,0;0,0;0,1];
    W = [W1,0;0,W2];
    C = expm([-A,gamma * W * gamma.';zeros(4),A.']);
    Ad = C(5:8,5:8).';
    Qd = Ad * C(1:4,5:8);
    x = zeros(4,length(xTrue));
    x(:,1) = x0;
    P = zeros(4,4,length(xTrue));
    P(:,:,1) = P0;
    K_plot = zeros(4,2,length(xTrue));
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
        K_plot(:,:,ind) = K;
        dx = K*dz;
        
        x(:,ind) = x(:,ind) + dx;
    end
end
