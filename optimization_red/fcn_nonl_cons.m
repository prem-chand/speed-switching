function [C,Ceq] = fcn_nonl_cons(ival,d,Fext,tF)

M_bez = 6;
x0 = ival(1:10)';
[alpha,theta_minus,theta_plus] = fcn_alpha_red(ival);
c = [1 1 0 0.5 0];
z0 = [c*x0(1:5)'; gamma0(x0)*x0(6:10)'];

options = odeset('Events',@(t,z) touchdown5_slope_red(t,z,alpha,theta_minus,theta_plus,d),'RelTol',1e-5,'AbsTol',1e-4);
[t,z,te,ze,ie] = ode45(@(t,z) fcn_dynamics_red(t,z,alpha,theta_minus,theta_plus),[0 20],z0,options);

if size(te,1) == 0
    Ceq = [];
    C = [];
else
    for tp=1:size(t)
        y(tp,:) = fcn_qu_to_q(z(tp,:),alpha,theta_minus,theta_plus);
        [fx,gx,gforce]=fcn_stance_dynamics5(y(tp,:));
        Fext_plot = [interp1(tF,Fext,t(tp));0];
        u(tp,:) = fcn_stance_controller5(y(tp,:),fx,gx,gforce,Fext_plot,alpha,theta_minus,theta_plus,M_bez);
        [h(tp,:),hdot(tp,:)] = fcn_liederivative(y(tp,:),alpha,theta_minus,theta_plus,M_bez);
        GRF(tp,:) = fcn_GRF(y(tp,:),u(tp,:));
        friction_cone(tp) = abs(GRF(tp,1)/GRF(tp,2));
    end
    q4 = y(:,4);
    q5 = y(:,5);
    ye = fcn_qu_to_q(ze,alpha,theta_minus,theta_plus);
    [xplus, Fimp] = impact_map5(ye);
%     gam0_plus = gamma0(xplus);
%     gam0_minus = gamma0(ye);
%     delta_zero = (gam0_plus*xplus(6:10))/(gam0_minus*ye(6:10));
%     zplus(1) = theta_plus;
%     zplus(2) = delta_zero*ye(2);
%     
%     zplus = [c*xplus(1:5); gamma0(xplus)*xplus(6:10)];
    friction_cone = [friction_cone abs(Fimp(1)/Fimp(2))];
    output_error = [h hdot];
    norm_output_error = max(sqrt(sum(output_error.^2,2)));
    step_length = fcn_position_swingfoot(ye);
%         initial_swing_foot=fcn_position_swingfoot(XX');
    average_speed=step_length(1)/te;

    
    C = [(max(max(abs(u)))-100)/100;
        -1*min(q4);
        -1*min(q5);
        (max(friction_cone) - 0.5)/0.5;
%         (norm_output_error - 1e-3)*1000;
        -1*(average_speed - 0.40)/0.40;
        (average_speed - 0.80)/0.80];

    Ceq(1) = norm(xplus(:)-x0(:));
end



end