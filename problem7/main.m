%% prob 4.
addpath("../scripts/functions")
addpath("../problem1/functions_p1")
close all

%% discretization
% continuous state matrix
A = [0,1,0,0;
    0,0,0,0;
    0,0,0,1;
    0,0,0,0];

% continuous Gamma Matrix (with Q_c = W)
gamma = [0,0;1,0;0,0;0,1];

% Gamma * W * Gamma^T
Q = gamma * gamma.';

% Van Loan's Discretization
C = expm([-A,Q;zeros(4), A.']);

Ad = C(5:8,5:8).';
Qd = Ad * C(1:4,5:8);

%% V and W
V0 = DEG2RAD^2;

R = eye(2) * V0;

W1 = 5.17E-4;
W2 = 2.67E-4;

Q = getQ(W1,W2);
%% true
[xTrue,t,z] = truetraj();
%% root locus
%{
dlqe:   x[n+1] = Ax + Bu + Gw
        y[n] = Cx + Du + v

[M,P,Z,E] = dlqe(A,G,C,Q,R)
where   Q = E[ww.']
        R = E[vv.']

thus we now have A,G,Q,R,
but C has to be calculated per every location
%}
drawlocus = @(R,Q,title_) drawlocusmain(xTrue,Ad,R,Q,title_);

drawlocus(R,Q,"R,Q")
drawlocus(R, getQ(100*W1,100*W2), "R,100Q")
drawlocus(100*R,getQ(W1,W2),"100R, Q")
drawlocus(100*R,getQ(100*W1,100*W2),"100R,100Q")

function drawlocusmain(xTrue,Ad,R,Q,title_)
    if nargin<5
        title_ = ""
    end
    figure
    f = waitbar(0,"root locus");
    for ind = 1:30:length(xTrue)
        waitbar(ind/length(xTrue),f,"calculating...");
        HTemp= Jacob_h(xTrue(:,ind));
        H = zeros(2,4);
        H(:,[1,3]) = HTemp;
        [~,~,~,E] = dlqe(Ad,eye(4),H,Q,R);
        
        for dim = 1:4
            rea = real(E(dim));
            ima = imag(E(dim));
            plot([rea,-rea],[ima,ima],'.k')
            drawnow
            hold on
        end
    end
    title(title_)
    close(f)
end
function Qd = getQ(W1,W2)
    % continuous state matrix
    A = [0,1,0,0;
        0,0,0,0;
        0,0,0,1;
        0,0,0,0];
    % continuous Gamma Matrix (with Q_c = W)
    gamma = [0,0; 1,0; 0,0; 0,1];
    
    C = expm([-A,gamma*[W1,0;0,W2]*gamma.';zeros(4),A.']);
    
    Ad = C(5:8,5:8).';
    Qd = Ad * C(1:4,5:8);
end