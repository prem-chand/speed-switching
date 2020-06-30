function clsloop=return_map5_cooperation(t,Z,stance_foot_zero,tF,betta,beta_c,beta_s,theta_minus,theta_plus,M_bez,N)
clsloop=zeros(10,1);
[fx,gx,gforce]=fcn_stance_dynamics5(Z);
[pL,pL_dot] = leader(t,Z);
Fext = compute_coop_force(Z,pL,pL_dot,stance_foot_zero);
u=fcn_stance_controller5_correction(Z,fx,gx,gforce,Fext,betta,beta_c,beta_s,theta_minus,theta_plus,M_bez,N);
% Fext
clsloop=fx+gx*u+gforce*Fext;
end
