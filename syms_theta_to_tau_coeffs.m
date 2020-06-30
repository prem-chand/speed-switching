% clear all;
syms theta theta_minus theta_plus tau
syms a0 a1 a2 a3 a4 a5 a6 a7 a8
p = a0 + a1*theta + a2*theta^2 + a3*theta^3 + a4*theta^4 + a5*theta^5 + a6*theta^6 + a7*theta^7 + a8*theta^8;
p_tau = simplify(subs(p,theta,tau*(theta_minus-theta_plus)+theta_plus));
p_tau = collect(p_tau,tau);
C = coeffs(p_tau,tau);
