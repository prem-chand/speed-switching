function eigenvalue_five_stair()

model_parameters;
load fixed_point_6.mat
x_fix     = fpStates(1:10);
alpha_red = fpStates(11:22);
alpha_s   = stance_param(x_fix,alpha_red);
p_stoe_initial = [0; 0];

% Integrator Constants
tfin = 5; relTol = 1e-7; absTol = 1e-7;
options   = odeset('Events',@guard_perturbation,...
    'RelTol',relTol,...
    'AbsTol',absTol,...
   'NormControl','off');
x = 1e-4; % perturbation

for i = 1:10
    x_fix = fpStates(1:10);
    %forward perturbation
    x_fix(i) = x_fix(i) + x;
    y0 = x_fix;
    x_s = [];
    [~,x_s,~,~,~]=...
        ode45(@stance_dynamics,[0,tfin],y0,options,alpha_s,p_stoe_initial);
    [M_s, ~] = size(x_s);
    pCOM = [];
    vCOM = [];
    for j = 1:M_s
    pCOM(j,:) = fcn_stance_pCOM(x_s(j,1:5));
    vCOM(j,:) = fcn_stance_vCOM(x_s(j,1:5), x_s(j,6:10));
    end
    x_s = [x_s(:,1:5) pCOM(:,1:2) x_s(:,6:10) vCOM(:,1:2)];
    % Impact Dynamics
    q_e_minus     = x_s(end,1:7);
    dq_e_minus    = x_s(end,8:14);
    [q_s,dq_s]    = stance_to_impact_transition(q_e_minus,dq_e_minus);
    y0_plus       = [q_s; dq_s];
    
    %backward perturbation
    x_fix(i) = x_fix(i)-2*x;
    y0 = x_fix;
    x_s = [];
    [~,x_s,~,~,~]=...
        ode45(@stance_dynamics,[0,tfin],y0,options,alpha_s,p_stoe_initial);
    [M_s, ~] = size(x_s);
    pCOM = [];
    vCOM = [];
    for j = 1:M_s
    pCOM(j,:) = fcn_stance_pCOM(x_s(j,1:5));
    vCOM(j,:) = fcn_stance_vCOM(x_s(j,1:5), x_s(j,6:10));
    end
    x_s = [x_s(:,1:5) pCOM(:,1:2) x_s(:,6:10) vCOM(:,1:2)];
    % Impact Dynamics
    q_e_minus     = x_s(end,1:7);
    dq_e_minus    = x_s(end,8:14);
    [q_s,dq_s]    = stance_to_impact_transition(q_e_minus,dq_e_minus);
    y0_minus      = [q_s; dq_s];
    
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
    
    