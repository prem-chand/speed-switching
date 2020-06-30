function check = compute_fp_trans_u_cone(speed_ctrl,ival,G,XX)
%COMPUTE_FP_SPEED_BETA Summary of this function goes here
%   Detailed explanation goes here
warning('OFF','ALL');
check = 1;
N = 8;
num_steps = 30;
M_bez=6;
te=0;
tplot=0;
tF=[0 500];
for i=1:num_steps
    clear theta_minus theta_plus betta event
%     [betta,beta_c,theta_minus,theta_plus]=fcn_alpha_red_correction(ival',XX(:,i),N);
    [betta,beta_c,theta_minus,theta_plus]=fcn_alpha_red_correction(ival',XX(:,1),N);
    beta_ctrl = pinv(G)*speed_ctrl;
    beta_s_p = fcn_beta_correction(XX(:,1),XX(:,i),beta_ctrl);
    for k=1:1:4
        beta_s_p_tau(k,:) = fcn_coeffs_theta_to_tau(beta_s_p(k,:),theta_plus,theta_minus);
    end
    Fext = [0 ; 0];
    options = odeset('Events',@(t,q)touchdown5_steps(t,q,0),'RelTol',1e-5,'AbsTol',1e-4);
    [t,y,te,ye,~] = ode45(@(t,Z) return_map5(t+tplot,Z,Fext,tF,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N),[0 20],XX(:,i),options);
    [XX(:,i+1) , ~] = impact_map5(ye);
    tplot=tplot+te;
    
    clear u friction_cone GRF
    for tp=1:size(t)
        [fx,gx,gforce]=fcn_stance_dynamics5(y(tp,:));
        u(tp,:)=fcn_stance_controller5_correction(y(tp,:)',fx,gx,gforce,Fext,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N);
        GRF(tp,:)=fcn_GRF(y(tp,:),u(tp,:));
        friction_cone(tp) = abs(GRF(tp,1)/GRF(tp,2));
    end
    fric_max = max(friction_cone);
    u_max = max(max(abs(u')));
    if (fric_max>0.81 || u_max>108.1)
        check = 0;
        break;
    end
end

end

