function cost=fixed_point_same_Z(para,initial_alpha_red,Fext,tF)
warning('off','all');

clear  theta_minus theta_plus betta XX event M_bez cost_d Fimp
% global theta_minus theta_plus betta event M_bez

XX=[initial_alpha_red(1:5); para(1:5)'];

%XX=initial_alpha_red(1:10);
M_bez=6; N=8;
[betta,beta_c,theta_minus,theta_plus]=fcn_alpha_red_correction([XX; initial_alpha_red(11:end)],XX,N);
beta_ctrl = 1e3*para(6:13)';
% beta_s_p = fcn_beta_correction(XX(:,1),XX(:,1),beta_ctrl);
beta_s_p = fcn_beta_correction(initial_alpha_red(1:10),XX(:,1),beta_ctrl);
for k=1:1:4
    beta_s_p_tau(k,:) = fcn_coeffs_theta_to_tau(beta_s_p(k,:),theta_plus,theta_minus);
end
options = odeset('Events',@(t,q) touchdown5_steps(t,q,0),'RelTol',1e-5,'AbsTol',1e-4);
[~,~,te,ye,~] = ode45(@(t,Z) return_map5(t,Z,Fext,tF,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N),[0 20],XX,options);
if  size(te,1)==0
    cost=20000;
else
    [xplus,~]=impact_map5(ye);
    cost=norm(xplus-XX);
end
end