function Dhc_Dq = fcn_Dhc_Dq(x,beta,theta_minus,theta_plus,N,tau_end)

c = [1 1 0 0.5 0];
tau = (c*x(1:5) - theta_plus)/(theta_minus - theta_plus);
if (tau>=tau_end || norm(norm(beta))<eps)
    Dhc_Dq = zeros(4,5);
else
    Dhc_Dq = zeros(4,5);
    Dhc_Dtau = zeros(4,1);
    for k=1:1:N
        Dhc_Dtau = Dhc_Dtau + (k*(tau^(k-1))*eye(4))*beta(:,k+1);
    end
    Dhc_Dq = (Dhc_Dtau*c)/(theta_minus - theta_plus);
end
end

