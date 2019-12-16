function [x_est,P,waypoint]=NewtonRhapson(x_est0,invM,Zmea,V)





    x_est_temp=x_est0;
    k=0;
    waypoint=x_est_temp;
    while(true)
        
        Z_est=h(x_est_temp);
%         J(k)=1/2*(x_est_temp-x_mean)'*invM*(x_est_temp-x_mean)+1/2*(Zmea(:,i)-Z_est)'*inv(V)*(Zmea(:,i)-Z_est);
        k=k+1;
        H_est=Jacob_h(x_est_temp);
        P=inv(invM+H_est'*inv(V)*H_est);
        GRD=invM*(x_est_temp-x_est0)-H_est'*inv(V)*(Zmea-Z_est);
        %norm(GRD)
        if(norm(GRD)<0.0000000001)
            %fprintf('Newton Rhapson iteration: %d',k);
            %k
            break;
        end
        x_est_temp=x_est_temp-P*GRD;
        waypoint=[waypoint,x_est_temp]; 
    end

    x_est=x_est_temp;
    



end