function clsloop=return_map5(t,Z,F,tF,betta,beta_c,beta_s,theta_minus,theta_plus,M_bez,N)
clsloop=zeros(10,1);
[fx,gx,gforce]=fcn_stance_dynamics5(Z);
Fext = [interp1(tF,F,t);0];
u=fcn_stance_controller5_correction(Z,fx,gx,gforce,Fext,betta,beta_c,beta_s,theta_minus,theta_plus,M_bez,N);
% Fext
clsloop=fx+gx*u+gforce*Fext;
end
