function dz = fcn_dynamics_red(t,z,alpha,theta_minus,theta_plus)

[x, kappa1] = fcn_qu_to_q(z,alpha,theta_minus,theta_plus);

% c = [1 1 0 0.5 0];
% H0 = [zeros(4,1) eye(4)];

% H = [H0; c];


% q = H\[hd; theta(1)];
% qdot = [dhd; gamma0(q)]\[zeros(4,1); 1];



dz = [kappa1*z(2);kappa2(x(1:5))];



















end