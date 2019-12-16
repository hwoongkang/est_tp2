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
x([2,4],1) = [50/3600;0];

x0 = x(:,1);

P0 = zeros(4);

P0([1,3],[1,3]) = P0pos;
P0([2,4],[2,4]) = 1E4 * eye(2);

P = zeros(4,4,length(xTrue));
P(:,:,1) = P0;

M = zeros(4,4,length(xTrue));

getXWithW = @(W) getX(Ad,Qd,W,R,xTrue,z,x0,P0);
getXWithW1AndW2 = @(W1,W2) getX2(Ad,Qd,W1,W2,R,xTrue,z,x0,P0);

x = getXWithW1AndW2(0.00025^2,0.0002^2);



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
        
        K = M(:,:,ind) * H.' * inv(H*M(:,:,ind) * H .' +R);
        P(:,:,ind) = (eye(4) - K*H) * M(:,:,ind);
        
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
