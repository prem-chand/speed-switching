function F = compute_coop_force(x,pL,pL_dot,stance_foot_zero)
%COMPUTE_COOP_FORCE Summary of this function goes here
%   Detailed explanation goes here
KL = 15;
NL = 10;
pE = fcn_position_head(x(1:5));
pE = pE'+[stance_foot_zero;0];
pE_dot = fcn_jacobian_force(x(1:5))*x(6:10);
% pL-pE
F = KL*(pL-pE) + NL*(pL_dot-pE_dot);

end

