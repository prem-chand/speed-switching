function speed = compute_speed_check_convergence(speed_ctrl,ival,G,XX)
%COMPUTE_FP_SPEED_BETA Summary of this function goes here
%   Detailed explanation goes here
warning('OFF','ALL');
N = 8;
num_steps = 30;
M_bez=6;
te=0;
tplot=0;
tF=[0 500];
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
    if length(ye~=0)
        [XX(:,i+1) , ~] = impact_map5(ye);
        tplot=tplot+te;
        step_length(i,:) = fcn_position_swingfoot(ye);
        average_speed(i) = step_length(i,1)/te;
    else
        speed = -1;
        break;
    end
end
speed = average_speed(end);

end

