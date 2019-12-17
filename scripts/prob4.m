%% prob 4.
addpath("functions")
addpath("../problem1/functions_p1")

%% discretization
% continuous state matrix
A = [0,1,0,0;
    0,0,0,0;
    0,0,0,1;
    0,0,0,0];

% continuous Gamma Matrix (with Q_c = W)
gamma = [0;1;0;1];
gamma2 = [0,0;1,0;0,0;0,1];

% Gamma * W * Gamma^T
Q = gamma * gamma.';
Q2 = gamma2 * gamma2.';
% Van Loan's Discretization
C = expm([-A,Q;zeros(4), A.']);
C2 = expm([-A,Q;zeros(4),A.']);

Ad = C(5:8,5:8).';
Ad2 = C2(5:8,5:8).';

Qd = Ad * C(1:4,5:8);
Qd2 = Ad2 * C2(1:4,5:8);

%% part (a)
[xTrue,t,z]= truetraj();

R = diag([DEG2RAD, DEG2RAD].^2);

[x0pos,P0pos]=NewtonRhapson([10;10],zeros(2),z(:,1),R);

% allocation
x = zeros(4,length(xTrue));
x([1,3],1) = x0pos;
x([2,4],1) = [0;0];

x0 = x(:,1);

P0 = zeros(4);

P0([1,3],[1,3]) = P0pos;
P0([2,4],[2,4]) = 1E10 * eye(2);

P = zeros(4,4,length(xTrue));
P(:,:,1) = P0;

M = zeros(4,4,length(xTrue));

%% function handles
% for (a) -> disturbability improved
getXWithW = @(W) getX(Ad,Qd + 3E-5*randn(4),W,R,xTrue,z,x0,P0);
% for (b)
getXWithW1AndW2 = @(W1,W2) getX2(A,Qd,W1,W2,R,xTrue,z,x0,P0);

% initial error var
err = 1E20;
WAns = 0;


f = waitbar(0,'W...');
for W = 8.9E-4:1E-6:9.1E-4%1E-5:1E-5:1E-3
    waitbar(W,f,sprintf("current best: %.3f^2",WAns));
    % get X
    x = getXWithW(W^2);
    % position mean square error
    poserr = x([1,3],:) - xTrue;
    errnow = trace(poserr.' * poserr) / length(x);
    
    % replace optimal
    if errnow<err
        err = errnow;
        WAns = W;
        bestXWithW = x;
    end
end
close(f)

% plotting
bestX = bestXWithW;
figure
plot(bestX(3,:),bestX(1,:),'.k')
hold on
plot(xTrue(2,:),xTrue(1,:),'--r')
title(sprintf("best result with W = %.2e",WAns))
axis equal


%% (b)
g = waitbar(0,'W1&W2...');

W1Ans = 0;
W2Ans = 0;
err2 = 1E20;

% best(W1, W2): (5.175, 2.67) * 1E-4
for W1 = 5.17*1E-4:1E-7:5.19*1E-4
    
    waitbar(W1*1000/2,g,sprintf("current best: %f,%f",W1Ans,W2Ans));
    
    for W2 = 2.67E-4 %2.6*1E-4:1E-6:2.8*1E-4
        
        % mean square error
        x2 = getXWithW1AndW2(W1^2,W2^2);
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

return


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

function x = getX2(A,Qd,W1,W2,R,xTrue,z,x0,P0)
    % van loan again
    gamma = [0,0;1,0;0,0;0,1];
    W = [W1,0;0,W2];
    C = expm([-A,gamma * W * gamma.';zeros(4),A.']);
    Ad = C(5:8,5:8).';
    Qd = Ad * C(1:4,5:8);
    % tuning parameter
    
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
