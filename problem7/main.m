%% prob 4.
addpath("../scripts/functions")
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

%% V and W
V0 = DEG2RAD^2;

R = eye(2) * V0;

W1 = 
Q = Qd * W0;
%% true
[xTrue,t,z] = truetraj();
%% root locus
%{
dlqe:   x[n+1] = Ax + Bu + Gw
        y[n] = Cx + Du + v

[M,P,Z,E] = dlqe(A,G,C,Q,R)
where   Q = E[ww.']
        R = E[vv.']

thus we now A,G,Q,R,
but C has to be calculated per every location
%}
for ind = 1%1:length(xTrue)
    HTemp= Jacob_h(xTrue(:,ind));
    H = zeros(2,4);
    H(:,[1,3]) = HTemp;
    [~,~,~,E] = dlqe(Ad,eye(4),H,Q,R);
end