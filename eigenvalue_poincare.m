function eigenvalue_poincare()

fix=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0.mat');
% ival=fix.fpStates
ival=fix.x'
x_fix     = ival(1:10);

clear  theta_minus theta_plus betta event XX M_bez
global theta_minus theta_plus betta event M_bez

[betta,theta_minus,theta_plus]=fcn_alpha_red(ival)
%event=x_fix(2);
M_bez=6;
tF=[0 50];
Fext=[0  0];

% Integrator Constants
tfin = 5; relTol = 1e-5; absTol = 1e-4;
options   = odeset('Events',@(t,q) touchdown5_steps(t,q,0),...
    'RelTol',relTol,...
    'AbsTol',absTol,...
   'NormControl','off');
x = 1e-4; % perturbation

for i = 1:10
    x_fix = ival(1:10);
    %forward perturbation
    x_fix(i) = x_fix(i) + x;
    y0 = x_fix;
    x_s = [];
    [~,~,~,ye,~]=...
        ode45(@(t,Z) return_map5_original(t,Z,Fext,tF,betta,theta_minus,theta_plus,M_bez),[0,tfin],y0,options);
    y0_plus = impact_map5(ye);
    
    %backward perturbation
    x_fix(i) = x_fix(i)-2*x;
    y0 = x_fix;
    x_s = [];
    [~,~,~,ye,~]=...
        ode45(@(t,Z) return_map5_original(t,Z,Fext,tF,betta,theta_minus,theta_plus,M_bez),[0,tfin],y0,options);
    y0_minus = impact_map5(ye);
    
    % A is the Jacobian matrix of poincare map at fixed point, compute the
    % i-th column of A
    A(1,i) = (y0_plus(1)-y0_minus(1))/(2*x);
    A(2,i) = (y0_plus(2)-y0_minus(2))/(2*x);
    A(3,i) = (y0_plus(3)-y0_minus(3))/(2*x);
    A(4,i) = (y0_plus(4)-y0_minus(4))/(2*x);
    A(5,i) = (y0_plus(5)-y0_minus(5))/(2*x);
    A(6,i) = (y0_plus(6)-y0_minus(6))/(2*x);
    A(7,i) = (y0_plus(7)-y0_minus(7))/(2*x);
    A(8,i) = (y0_plus(8)-y0_minus(8))/(2*x);
    A(9,i) = (y0_plus(9)-y0_minus(9))/(2*x);
    A(10,i) = (y0_plus(10)-y0_minus(10))/(2*x);    
end
A;
d=eig(A)
for i=1:1:10
    abs(d(i))
end
    
    