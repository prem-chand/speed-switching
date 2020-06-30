function [WF,diff] = compute_wf_bound(speed_ctrl,ival,G,XX)
%COMPUTE_FP_SPEED_BETA Summary of this function goes here
%   Detailed explanation goes here
warning('OFF','ALL');
c = [1 1 0 0.5 0];
N = 8;
num_steps = 1;
M_bez=6;
te=0;
tF=[0 500];
tplot = 0;
F = 0.95;
for i=1:num_steps
    clear theta_minus theta_plus betta event
    
    [betta,beta_c,theta_minus,theta_plus]=fcn_alpha_red_correction(ival',XX(:,i),N);
    beta_ctrl = pinv(G)*speed_ctrl;
    beta_s_p = fcn_beta_correction(XX(:,1),XX(:,i),beta_ctrl);
    for k=1:1:4
        beta_s_p_tau(k,:) = fcn_coeffs_theta_to_tau(beta_s_p(k,:),theta_plus,theta_minus);
    end
    Fext = [0 ; 0];
    options = odeset('Events',@(t,q)touchdown5_steps(t,q,0),'RelTol',1e-5,'AbsTol',1e-4);
    [t,y,te,ye,~] = ode45(@(t,Z) return_map5(t+tplot,Z,Fext,tF,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N),[0 20],XX(:,i),options);
    [XX(:,i+1) , ~] = impact_map5(ye);
    gam0_plus=gamma0(XX(:,i+1));
    gam0_minus=gamma0(ye);
    delta_zero=(gam0_plus*XX(6:10,i+1))/(gam0_minus*ye(6:10)');
    
    for tp=1:size(t)
        theta(tp)=c*y(tp,1:5)';
        s(tp)=(theta(tp)-theta_plus)/(theta_minus-theta_plus);
        k1(tp)=c*([fcn_Dh_Dq(s(tp),betta,(theta_plus-theta_minus))-fcn_Dhc_Dq(y(tp,:)',beta_s_p_tau,theta_minus,theta_plus,7,0.9)-fcn_Dhc_Dq(y(tp,:)',beta_c,theta_minus,theta_plus,N,0.5) ; gamma0(y(tp,1:5))]\[zeros(4,1) ; 1]);
        k2(tp)=kappa2(y(tp,1:5));
        ksi1(tp)=theta(tp);
        ksi2(tp)=gamma0(y(tp,1:5)')*y(tp,6:10)';
        k3(tp,:)=kappa3(y(tp,1:5));
        w(tp)=norm(k3(tp,:))/k1(tp);
    end
    
    Vzero=trapz(theta,-k2./k1);
    WF=-F*trapz(theta,w); 
    for j=2:1:length(t)
        Vzero_ksi1(j-1)=trapz(theta(1:j),-k2(1:j)./k1(1:j));
        WF_ksi1(j-1)=-F*trapz(theta(1:j),w(1:j));
    end
    M=max(Vzero_ksi1-WF_ksi1);
    K=max(Vzero_ksi1);
%     diff=WF-Vzero-((1-delta_zero^2)/delta_zero^2)*M;
%     
%     zeta_minus_lb = K/(delta_zero^2);
    zeta_minus_fp = 0.5*(gam0_minus*ye(6:10)')^2;
%     zeta_minus_fp
%     -Vzero/(1-delta_zero^2)
%     diff = zeta_minus_fp + WF/(1-delta_zero^2)-M/delta_zero^2
    diff=WF-Vzero-((1-delta_zero^2)/delta_zero^2)*M;
%     if i>1
%         norm(XX(:,i+1)-XX(:,i+1))
%     end
    tplot=tplot+te;
end

end

