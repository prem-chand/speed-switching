function u=fcn_stance_controller5(Y,f,g,gforce,Fext,betta,theta_minus,theta_plus,M_bez)  %Y=state,f=f(x),g=g(x),gforce=gf(x)
% global betta M_bez

% if M_bez==6 
    [h,Lfh,Lfhx]=fcn_liederivative(Y,betta,theta_minus,theta_plus,M_bez);
% elseif M_bez==5
%     [h,Lfh,Lfhx]=fcn_liederivative_m5(Y,betta);
% end

Lffh=Lfhx*f;
LgLfh=Lfhx*g;  %4*4
LgforceLfh=Lfhx*gforce; %4*2
%PD controller design
% alpha=0.9;
% eps=0.05;
% fai=@(x1,x2) x1+(0.5-alpha)*sign(x2)*abs(x2)^(2-alpha);
% v=@(x1,x2) -sign(x2)*abs(x2)^alpha-sign(fai(x1,x2))*abs(fai(x1,x2))^(alpha/(2-alpha));
% sai=[1/eps^2*v(h(1),eps*Lfh(1));1/eps^2*v(h(2),eps*Lfh(2));1/eps^2*v(h(3),eps*Lfh(3));1/eps^2*v(h(4),eps*Lfh(4))];

% sai=0;
%control law

eta = [h ; Lfh];
LgLfh_inv = inv(LgLfh);
epsilon=0.1;
Kp=16*eye(4);
Kd=8*eye(4);

%% PD Controller
% sai=-Kp*h*(1/epsilon)^2-Kd*Lfh*(1/epsilon);
% u=LgLfh\(sai-Lffh-LgforceLfh*Fext);

%% QP Controller
A = [zeros(4,4) eye(4,4) ; -Kp -Kd];
Q = eye(8);
P = lyap(A',Q);
I_eps = [(1/epsilon)*eye(4) zeros(4) ; zeros(4) eye(4)];
P_eps = I_eps*P*I_eps;
V_eps = eta'*P_eps*eta;
[~,eig_Q] = eig(Q);
[~,eig_P] = eig(P);
c3 = min(eig(Q))/max(eig(P));
F = [zeros(4) eye(4) ; zeros(4) zeros(4)];
G = [zeros(4) ; eye(4)];
psi_0 = eta'*(F'*P_eps + P_eps*F)*eta + (c3/epsilon)*V_eps;
psi_1 = (2*eta'*P_eps*G);
u_min = -80;
u_max = 80;
u_star = -LgLfh_inv*Lffh;
p1 = 1000;
H = [2*eye(4) zeros(4,1) ; zeros(1,4) 2*p1]; % matrix of the quadratic cost, p1 is the penalty
f_lin = zeros(1,5)';
A_ineq = [psi_1 -1 ; -LgLfh_inv zeros(4,1) ; LgLfh_inv zeros(4,1)];
B_ineq = [-psi_0 ; u_star-u_min ; u_max-u_star];

options = optimoptions('quadprog','Display','off');
optim_sol = quadprog(H,f_lin,A_ineq,B_ineq,[],[],[],[],[],options);
mu = optim_sol(1:4,1);
u = u_star + (LgLfh_inv)*mu;


% Fext
% u=LgLfh\(sai-Lffh-LgforceLfh*Fext);
%lgf_Force=LgforceLfh*Fext;
%ydotdot=Lffh+LgLfh*u;
end