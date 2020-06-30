function [value,isterminal,direction] = touchdown5(t,q,angle)
% Locate the time when height passes through zero in a 
% decreasing direction and stop integration.
% global event
% global theta_minus
% c=[1 1 0 0.5 0];
% value=theta_minus-c*q(1:5);
% value=q(2)-XX(3,1);
% value=q(3)-XX(2,1)

% value=q(3)-event;
% angle=-2; 
%Note - Change angle in compute_nonlcon also
Psfoot=fcn_position_swingfoot(q);
Psfoot_vel=fcn_swing_foot_velocity_vert(q);
% value=Psfoot(2);
vel=norm(Psfoot_vel)*cos(atan2(Psfoot_vel(1),Psfoot_vel(2))+angle*pi/180);
value=1;
if (Psfoot(2)/Psfoot(1)<tand(angle)) && (vel<-0.02) && Psfoot(1)>0
   value=-1;
end
isterminal = 1;   % Stop the integration
direction = -1;   % Negative direction only
end