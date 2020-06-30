function clsloop=return_map5_original(t,Z,Fext,tF,betta,theta_minus,theta_plus,M_bez)
clsloop=zeros(10,1);
[fx,gx,gforce]=fcn_stance_dynamics5(Z);
Fext = [interp1(tF,Fext,t);0];
u=fcn_stance_controller5(Z,fx,gx,gforce,Fext,betta,theta_minus,theta_plus,M_bez);
% Fext
clsloop=fx+gx*u+gforce*Fext;
end
