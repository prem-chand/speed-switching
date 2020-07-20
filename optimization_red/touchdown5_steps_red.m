function [value,isterminal,direction] = touchdown5_steps_red(t,qu,alpha,theta_minus,theta_plus,d)

x = fcn_qu_to_q(qu,alpha,theta_minus,theta_plus);

Psfoot=fcn_position_swingfoot(x(1:5));


value = Psfoot(2) - d;


isterminal = 1;   % Stop the integration
direction = -1;   % Negative direction only
end