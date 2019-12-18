clc;
clear;
close all;
%% prob 5.

addpath("problem1/functions_p1")

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

%%
[xTrue,t,z,vTrue]= truetraj();
stateTrue=[xTrue(1,:);vTrue(1,:);xTrue(2,:);vTrue(2,:)];

R = diag([DEG2RAD, DEG2RAD].^2);

[x0pos,P0pos]=NewtonRhapson([10;10],zeros(2),z(:,1),R);

%% initial change

x0pos_10 = xTrue(:,1) + 100*(x0pos - xTrue(:,1)); %%%%%%%%%%%%

% allocation

x=zeros(4,1);

x([1,3],1) = x0pos_10;
x([2,4],1) = [50/3600;0];

x0 = x(:,1);

P0 = zeros(4);

P0pos_10 = P0pos * (100)^2;  %%%%%%%%%%%%

P0([1,3],[1,3]) = P0pos_10;
P0([2,4],[2,4]) = 1E4 * eye(2);

M = zeros(4,4,length(xTrue));

getXKPrWithW1AndW2 = @(W1,W2) getXKPr2(Ad,Qd,W1,W2,R,xTrue,z,x0,P0);


%% part(b)



g = waitbar(0,'W1&W2...');

W1Ans = 0;
W2Ans = 0;
err2 = 1E20;

% best(W1, W2): (5.18, 2.67) * 1E-4
for W1 = 0.0228^2 %[0.001^2, 0.0228^2 1]%[0.00001^2, 0.0228^2, 0.05^2, 0.1^2   ]%5.17*1E-4:1E-6:5.19*1E-4
    
    waitbar(W1*1000/2,g,sprintf("current best: %f,%f",W1Ans,W2Ans));
    
    for W2 = 0.0163^2 %[0.001^2, 0.0163^2 1]%[0.001^2, 0.0163^2, 0.05^2, 0.1^2]  % 2.6*1E-4:1E-5:2.8*1E-4
        
        [x2,K_list,P,r] = getXKPrWithW1AndW2(W1^2,W2^2);
        poserr2 = x2([1,3],:) - xTrue;
        stateErr= abs(x2-stateTrue);
        errnow = trace(poserr2.' * poserr2)/length(x2);
        if errnow<err2
            bextXWithW1AndW2 = x2;
            err2 = errnow;
            W1Ans = W1;
            W2Ans = W2;
        end
        
        figure;sgtitle(sprintf("5(b) Kalman Gain W1=%.5f^2 W2=%.5f^2",sqrt(W1),sqrt(W1)));
        for i=1:4
            for j=1:2
                subplot(4,2,(i-1)*2+j);plot(t(2:end),K_list(:,i,j));title(sprintf("K(%d,%d)",i,j));
            end
        end
      
        figure(1);hold on; title(sprintf("6 Trajectory W1=%.5f^2 W2=%.5f^2",sqrt(W1),sqrt(W2)));
        f1=plot(x2(3,:),x2(1,:),'k'); hold on; f2=plot(xTrue(2,:),xTrue(1,:),'r','LInewidth',3);f3=plot(x2(3,1),x2(1,1),'k^');
        legend([f1,f2,f3],'est','true','start');
        axis equal;
       
        figure; sgtitle(sprintf("6 Err/Cov W1=%.5f^2 W2=%.5f^2",sqrt(W1),sqrt(W2)));
        for i=1:4 
            subplot(4,1,i);plot(3*sqrt(P(1:end,i,i)),'r'); hold on; plot(stateErr(i,1:end),'k');
%         if(i==2 || i==4) ylim([0 2]); end
            if( i==1) legend('3\sigma','error'); end
        end
        
        figure; sgtitle(sprintf("6 Residual W1=%.5f^2 W2=%.5f^2",sqrt(W1),sqrt(W2)));
        for i=1:2
            subplot(2,1,i); hold on; plot(abs(r(i,:)),'k');
        end
        
    end
end

fprintf("4(b) best Answer W1=%f W2=%f \n",W1Ans,W2Ans);

close(g)
%% Functions

function [x,K_list,P,r] = getXKPr2(Ad,Qd,W1,W2,R,xTrue,z,x0,P0)
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
    P = zeros(length(xTrue),4,4);
    r= zeros(2,length(xTrue)-1);
    P(1,:,:) = P0;
    
%     Ptemp=reshape(P(1,:,:),4,4);
%     [evector,eval]=eigs(Ptemp([1,3],[1,3]));
%     eli=elipse(evector,3.*sqrt(eval),x([1,3],1));
%     figure(1); hold on;
%     plot(eli(2,:),eli(1,:),'b');
    for ind = 2:length(xTrue)
        % time update
        x(:,ind) = Ad * x(:,ind-1);
        Ptemp=reshape(P(ind-1,:,:),4,4);
        M = Ad * Ptemp * Ad.' + Qd;
        
        % measurement update
        dz = z(:,ind) - h(x([1,3],ind));
        r(:,ind-1)=dz;
        Htemp = Jacob_h(x([1,3],ind));
        H = zeros(2,4);
        H(:,[1,3]) = Htemp;
        
        K = M * H.' * inv(H*M * H .' +R);
        P(ind,:,:) = (eye(4) - K*H) * M;
        
             
        dx = K*dz;
        
        K_list(ind-1,:,:)=K;
        x(:,ind) = x(:,ind) + dx;
        
%         Ptemp=reshape(P(ind,:,:),4,4);
%         [evector,eval]=eigs(Ptemp([1,3],[1,3]));
%         eli=elipse(evector,3.*sqrt(eval),x([1,3],ind));
%         figure(1); hold on;
%         plot(eli(2,:),eli(1,:),'b');
    end
end