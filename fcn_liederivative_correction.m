function [hc,Lfhc,DLfhc_Dx] = fcn_liederivative_correction(x,beta,theta_minus,theta_plus,N,tau_end)

%   global theta_minus theta_plus;

c = [1 1 0 0.5 0];
tau = (c*x(1:5) - theta_plus)/(theta_minus - theta_plus);
tau_dot = (c*x(6:10))/(theta_minus - theta_plus);
if (tau>=tau_end || norm(norm(beta))<eps)
    hc = zeros(4,1);
    Lfhc = zeros(4,1);
    DLfhc_Dx = zeros(4,10);
else
    hc = zeros(4,1);
    Lfhc = zeros(4,1);
    DLfhc_Dx = zeros(4,10);
    for k=0:1:N
        hc = hc + ((tau^k)*eye(4))*beta(:,k+1);
    end
    
    Dhc_Dtau = zeros(4,1);
    for k=1:1:N
        Dhc_Dtau = Dhc_Dtau + (k*(tau^(k-1))*eye(4))*beta(:,k+1);
    end
    Lfhc = Dhc_Dtau*tau_dot;

    D2hc_Dtau2 = zeros(4,1);
    for k=2:1:N
        D2hc_Dtau2 = D2hc_Dtau2 + (k*(k-1)*(tau^(k-2))*eye(4))*beta(:,k+1);
    end
    DLfhc_Dx = [D2hc_Dtau2*tau_dot*c Dhc_Dtau*c]/(theta_minus-theta_plus);
end
end

