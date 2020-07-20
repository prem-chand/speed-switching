clear;
clc; 

tF = [0 20];
Fext = [0 0];
d = 2.5;

fix = open('optimization_red/fp1-2p5.mat');
% fix=open('fixedpointforfivelink/Switching_Control/fixed_point_d=1cm.mat');
ival = fix.x';
% ival(11:end) = ival(11:end)*1.0001;
[alpha, theta_minus, theta_plus] = fcn_alpha_red(ival);

% x0 = ival(1:10)';
% c = [1 1 0 0.5 0];
% z0 = [c*x0(1:5)'; gamma0(x0)*x0(6:10)'];

% LB = [ival(1:10) - 0.1; ival(11:end) - 0.1];
% UB = [ival(1:10) + 0.1; ival(11:end) + 0.1];

LB = [ival(1:5) - 0.00001; ival(6:10) - 0.1; ival(11:end) - 0.1];
UB = [ival(1:5) + 0.00001; ival(6:10) + 0.1; ival(11:end) + 0.1];

options = optimset('Display','iter','UseParallel',false,'TolFun',1e-4,'TolX',1e-10,'MaxFunEvals',30000);
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(z) cost_fn_red(z,d),ival,[],[],[],[],LB,UB,@(z) fcn_nonl_cons(z,d,Fext,tF),options)

cost_verify = cost_fn_red(x,d)
fcn_nonl_cons(x,d,Fext,tF)