clear;
clc; 

M=6;lf=0.4;lt=0.4;
clear  theta_minus theta_plus betta XX event
tF=[0 20];
Fext=[0 0];
d = 1;

% fix=open('fixedpointforfivelink/Switching_Control/fixed_point_d=1cm.mat');
fix = open('optimization_red/fp1-1.mat');
% ival1=fix.fpStates;
ival1=fix.x';

LB=[ival1(1:10)-0.1 ; ival1(11:end)-0.1];
UB=[ival1(1:10)+0.1 ; ival1(11:end)+0.1];

options = optimset('Display','iter','UseParallel',true','TolFun',1e-4,'TolX',1e-8,'MaxFunEvals',30000);
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(para) fixed_point2(para,Fext,tF,d),ival1',[],[],[],[],LB,UB,@(para) compute_nonlcon(para,Fext,tF,d),options)

cost_verify = fixed_point2(x,Fext,tF,d)
compute_nonlcon(x,Fext,tF,d)

% save('fixedpointforfivelink/Switching_Control/fixed_point_F=0_slow.mat','x');
