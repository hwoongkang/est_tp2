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
x([2,4],1) = [0;0];

x0 = x(:,1);

P0 = zeros(4);

P0([1,3],[1,3]) = P0pos;
P0([2,4],[2,4]) = 1E4 * eye(2);

P = zeros(4,4,length(xTrue));
P(:,:,1) = P0;

M = zeros(4,4,length(xTrue));

getXWithW = @(W) getX(Ad,Qd,W,R,xTrue,z,x0,P0);
getXWithW1AndW2 = @(W1,W2) getX2(Ad,Qd,W1,W2,R,xTrue,z,x0,P0);

err = 1E20;
WAns = 0;
f = waitbar(0,'W...');



for W = 8.9E-4:1E-6:9.1E-4%1E-5:1E-5:1E-3
    waitbar(W,f,sprintf("current best: %.3f^2",WAns));
    x = getXWithW(W^2);
    poserr = x([1,3],:) - xTrue;
    errnow = trace(poserr.' * poserr) / length(x);
    if errnow<err
        err = errnow;
        WAns = W;
        bestXWithW = x;
    end
end
WAns
close(f)

x = getXWithW1AndW2(0.00025^2,0.0002^2);

g = waitbar(0,'W1&W2...');

W1Ans = 0;
W2Ans = 0;
err2 = 1E20;

% best(W1, W2): (5.18, 2.67) * 1E-4
for W1 = 5.17*1E-4:1E-7:5.19*1E-4
    
    waitbar(W1*1000/2,g,sprintf("current best: %f,%f",W1Ans,W2Ans));
    
    for W2 = 2.6*1E-4:1E-6:2.8*1E-4
        
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

return
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

function x = getX2(Ad,Qd,W1,W2,R,xTrue,z,x0,P0)
    % tuning parameter
    Qd = zeros(4);
    Qd(2,2) = W1;
    Qd(4,4) = W2;
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
