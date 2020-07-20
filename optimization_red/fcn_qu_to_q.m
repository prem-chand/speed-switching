function [x,kappa1] = fcn_qu_to_q(z,alpha,theta_minus,theta_plus)

theta = z(1);

[hd,dhd] = desired_outputs(theta,alpha,theta_minus,theta_plus);

c = [1 1 0 0.5 0];
H0 = [zeros(4,1) eye(4)];

H = [H0; c];

dh = H0 - dhd'*c;
q = H\[hd'; theta(1)];
a = [dh; gamma0(q)]\[zeros(4,1); 1];
qdot = a*z(2);
kappa1 = c*a;
x = [q; qdot];



end