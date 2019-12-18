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
drawgif = @(R,Q,title_,todraw) drawlocusgif(xTrue,Ad,R,Q,title_,todraw);
drawmain = @(R,Q,title_) drawlocusmain(xTrue,Ad,R,Q,title_);
% % % drawgif(R,Q,"R,Q",true)
% % % drawgif(R, getQ(100*W1,100*W2), "R,100Q",1)
% % % drawgif(100*R,getQ(W1,W2),"100R, Q",1)
% % % drawgif(100*R,getQ(100*W1,100*W2),"100R,100Q",1)

drawmain(R,Q,"R, Q")
drawmain(R, getQ(100*W1,100*W2), "R, 100Q")
drawmain(100*R,getQ(W1,W2),"100R, Q")
drawmain(100*R,getQ(100*W1,100*W2),"100R, 100Q")

function drawlocusgif(xTrue,Ad,R,Q,title_,todraw)
    filename = sprintf("%s.gif",strrep(title_,",","_"));
    if nargin<6
        todraw = false;
    end
    if nargin<5
        title_ = "";
    end
    
    % bigger characters
    fig = figure('DefaultAxesFontSize',20,'DefaultLineLineWidth',1);
    % set size
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.2, 0.7, 0.5]);
    
    subplot(1,2,1)
    title("Trajectory")
    plot(xTrue(2,:),xTrue(1,:),'--r');
    subplot(1,2,2)
    title("Root locus")
    sgtitle(title_)
    
    % waitbar
    f = waitbar(0,"root locus");
    
    % save previous eigenvalues
    savE = zeros(2,8*length(xTrue));
    % for convenience
    count =-8;
    % initial work needed
    framecount=1;
    for ind = 2:60:length(xTrue)
        % 4 eigenvalues, symmetric -> 8 ev's for each iteration
        count = count+8;
        waitbar(ind/length(xTrue),f,"calculating...");
        
        % draw the whole traj and current location
        subplot(1,2,1)
        hold off
        plot(xTrue(2,:),xTrue(1,:),'--r')
        hold on
        plot(xTrue(2,ind),xTrue(1,ind),'sk','MarkerSize',12,'MarkerFaceColor','k');
        drawnow
        title("Trajectory")
        
        % calculating the eigenvalues
        HTemp= Jacob_h(xTrue(:,ind));
        H = zeros(2,4);
        H(:,[1,3]) = HTemp;
        [~,~,~,E] = dlqe(Ad,eye(4),H,Q,R);
        
        % put eigenvalues to the history
        for dim = 1:4
            rea = real(E(dim));
            ima = imag(E(dim));
            savE(1,count+2*dim) = rea;
            savE(1,count+2*dim-1) = -rea;
            savE(2,count+2*dim) = ima;
            savE(2,count+2*dim-1) = ima;
            %plot([rea,-rea],[ima,ima],'.k')
            %drawnow
            %hold on
        end
        
        % draw the history of eigenvalues and annotate current ones
        subplot(1,2,2)
        hold off
        plot(savE(1,1:count+8),savE(2,1:count+8),'.k');
        hold on
        plot(savE(1,count+1:count+8),savE(2,count+1:count+8),'or','MarkerSize',7,'MarkerFaceColor','r')
        title("Symmetric Root Locus")
        drawnow
        
        if todraw
            % writing gif
            frame = getframe(fig);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            
            % initial work, 25fps
            if framecount==1
                imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0.04);
            % after the first frame, 25fps
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.04);
            end
            framecount = framecount+1;
        end
        
        % !! no saving process needed!
    end
    
    close(f)
end
function drawlocusmain(xTrue,Ad,R,Q,title_)
    if nargin<5
        title_ = ""
    end
    fig = figure('DefaultAxesFontSize',32,'DefaultLineLineWidth',1);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.2, 0.5, 0.6]);
    f = waitbar(0,"root locus");
    for ind = 1:60:length(xTrue)
        waitbar(ind/length(xTrue),f,"calculating...");
        HTemp= Jacob_h(xTrue(:,ind));
        H = zeros(2,4);
        H(:,[1,3]) = HTemp;
        [~,~,~,E] = dlqe(Ad,eye(4),H,Q,R);
        
        for dim = 1:4
            rea = real(E(dim));
            ima = imag(E(dim));
            plot([rea,-rea],[ima,ima],'ok','MarkerSize',5,'MarkerFaceColor','k')
            drawnow
            hold on
        end
    end
    title(title_)
    grid on
    saveas(gcf,sprintf("../figures/%s.jpg",strrep(title_,",","_")))
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