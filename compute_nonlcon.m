function [C,Ceq]=compute_nonlcon(initial_alpha_red,Fext,tF)

clear  theta_minus theta_plus betta XX event output_error
M_bez=6;
d = 0;
angle = 0;
XX=initial_alpha_red(1:10)';    % for full 26 parameters optimization

[betta,theta_minus,theta_plus]=fcn_alpha_red(initial_alpha_red');
count1=0; count2=0;

for i = 1:1:length(d)
    options = odeset('Events',@(t,q) touchdown5_steps(t,q,d(i)),'RelTol',1e-5,'AbsTol',1e-4);
    [t,y,te,ye,ie] = ode45(@(t,Z) return_map5_original(t,Z,Fext,tF,betta,theta_minus,theta_plus,M_bez),[0 20],XX,options);
    
    if size(te,1)==0
        Ceq=[];
        C=[];
    else
        for tp=1:size(t)
            [fx,gx,gforce]=fcn_stance_dynamics5(y(tp,:));
            Fext_plot = [interp1(tF,Fext,t(tp));0];
            u(tp,:) = fcn_stance_controller5(y(tp,:),fx,gx,gforce,Fext_plot,betta,theta_minus,theta_plus,M_bez);
            [h(tp,:),hdot(tp,:)] = fcn_liederivative(y(tp,:),betta,theta_minus,theta_plus,M_bez);
            GRF(tp,:) = fcn_GRF(y(tp,:),u(tp,:));
            friction_cone(tp) = abs(GRF(tp,1)/GRF(tp,2));
        end
        q4 = y(:,4);
        q5 = y(:,5);
        [fp , Fimp] = impact_map5(ye);
        friction_cone = [friction_cone abs(Fimp(1)/Fimp(2))];
        output_error=[h hdot];
        norm_output_error = max(sqrt(sum(output_error.^2,2)));
        step_length=fcn_position_swingfoot(ye);
%         initial_swing_foot=fcn_position_swingfoot(XX');
        average_speed=step_length(1)/te;
%         C=[(max(max(abs(u)))-200)/200 ; abs((norm_output_error-1e-3))*10000 ; (max(friction_cone)-0.7)/0.7 ; -1*min(q4) ; -1*min(q5);...
%         (average_speed+0.10)/0.10];
        C=[(max(max(abs(u)))-100)/100 ; -1*min(q4) ; -1*min(q5); (max(friction_cone)-0.5)/0.5 ; (norm_output_error-1e-3)*1000; ...
            -1*(average_speed-0.40)/0.40 ; (average_speed-0.50)/0.50];
        Ceq=norm(fp-XX);
    end
end
% maximum allowed torque set to 60
%   size(u)
%   C=max(max(abs(u),[],2))-100;
% speed range set
%C=[0.35-average_speed;average_speed-0.6]

%Ceq=average_speed-1.1
end