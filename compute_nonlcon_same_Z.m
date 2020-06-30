function [C,Ceq]=compute_nonlcon_same_Z(para,initial_alpha_red,Fext,tF)

clear  theta_minus theta_plus betta XX event output_error
M_bez=6;
XX=[initial_alpha_red(1:5); para(1:5)'];

%XX=initial_alpha_red(1:10);
M_bez=6; N=8;
[betta,beta_c,theta_minus,theta_plus]=fcn_alpha_red_correction(initial_alpha_red,XX,N);
beta_ctrl = 1e3*para(6:13)';
beta_s_p = fcn_beta_correction(XX(:,1),XX(:,1),beta_ctrl);
for k=1:1:4
    beta_s_p_tau(k,:) = fcn_coeffs_theta_to_tau(beta_s_p(k,:),theta_plus,theta_minus);
end
options = odeset('Events',@(t,q) touchdown5_steps(t,q,0),'RelTol',1e-5,'AbsTol',1e-4);
[t,y,te,ye,ie] = ode45(@(t,Z) return_map5(t,Z,Fext,tF,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N),[0 20],XX,options);

if size(te,1)==0
    Ceq=[];
    C=[];
else
    for tp=1:size(t)
        [fx,gx,gforce]=fcn_stance_dynamics5(y(tp,:));
        Fext_plot = [interp1(tF,Fext,t(tp));0];
        u(tp,:)=fcn_stance_controller5_correction(y(tp,:)',fx,gx,gforce,Fext,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N);
        [h(tp,:),hdot(tp,:)] = fcn_liederivative(y(tp,:),betta,theta_minus,theta_plus,M_bez);
        [hc(tp,:),hcdot(tp,:),~]=fcn_liederivative_correction(y(tp,:)',beta_c,theta_minus,theta_plus,N,0.5);
        [hs(tp,:),hsdot(tp,:),~]=fcn_liederivative_correction(y(tp,:)',beta_s_p_tau,theta_minus,theta_plus,7,0.9);
        GRF(tp,:) = fcn_GRF(y(tp,:),u(tp,:));
        friction_cone(tp) = abs(GRF(tp,1)/GRF(tp,2));
    end
    q4 = y(:,4);
    q5 = y(:,5);
    [fp , Fimp] = impact_map5(ye);
    friction_cone = [friction_cone abs(Fimp(1)/Fimp(2))];
    output_error=[h-hc-hs hdot-hcdot-hsdot];
    norm_output_error = max(sqrt(sum(output_error.^2,2)));
    step_length=fcn_position_swingfoot(ye);
    %         initial_swing_foot=fcn_position_swingfoot(XX');
    average_speed=step_length(1)/te;
    %         C=[(max(max(abs(u)))-200)/200 ; abs((norm_output_error-1e-3))*10000 ; (max(friction_cone)-0.7)/0.7 ; -1*min(q4) ; -1*min(q5);...
    %         (average_speed+0.10)/0.10];
    C=[(max(max(abs(u)))-200)/200 ; -1*min(q4) ; -1*min(q5); (max(friction_cone)-0.5)/0.5 ; (norm_output_error-1e-3)*1000];
    
%     C=[(max(max(abs(u)))-100)/100 ; -1*min(q4) ; -1*min(q5); (max(friction_cone)-0.5)/0.5 ; (norm_output_error-1e-3)*1000; ...
%         -1*(average_speed-1.00)/1.00];
    Ceq=norm(fp-XX);
end
end
% maximum allowed torque set to 60
%   size(u)
%   C=max(max(abs(u),[],2))-100;
% speed range set
%C=[0.35-average_speed;average_speed-0.6]

%Ceq=average_speed-1.1