%% prob 4.
addpath("functions")
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

% filtering
for ind = 1:length(xTrue)
end