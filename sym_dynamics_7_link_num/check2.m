clear all;
xl=rand([10 1]);
x=[xl(1:5) ; 1-xl(1) ; 1 ; xl(6:10) ; -xl(6) ; 0];
[Gl,Dl,Cl]=fcn_dynamics_matrices2_const(xl);
[G,D,C]=Copy_of_fcn_dynamics_matrices2(x);
D_l=[D(1:5,1)-D(1:5,6) D(1:5,2:5)];
C_l=[C(1:5,1)-C(1:5,6) C(1:5,2:5)];
G_l=G(1:5);
ddq=rand([5 1]);
LHSl=Dl*ddq+Cl*xl(6:10)+Gl;
LHS_l=D_l*ddq+C_l*xl(6:10)+G_l;
LHSl-LHS_l