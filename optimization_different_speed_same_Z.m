clear;
clc; 

M=6;lf=0.4;lt=0.4;
clear  theta_minus theta_plus betta XX event
tF=[0 20];
Fext=[5;0];

% fix=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0.mat');
fix=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0_slow.mat');
% ival1=fix.fpStates;
fp=fix.x';

fix2 = open('fixedpointforfivelink/Switching_Control/fixed_point_F=5_same_Z.mat');

% beta_ctrl_guess = 1e-3*1*1.0e+06*[0.5783   -1.0944   -2.4011    4.5318   -0.6196    1.1696    1.7663   -3.3349]';
% beta_ctrl_guess = ones(8,1);
beta_ctrl_guess = zeros(8,1);

ival1 = [fp(6:10); beta_ctrl_guess]*1
% ival1 = [fix2.x_save(6:10) ; fix2.x_save(27:end)]*1.1;

LB=[ival1(1:5)-0.1 ; ival1(6:end)-1000] ;
UB=[ival1(1:5)+0.1 ; ival1(6:end)+1000];

% LB=[ival1(1:10)-0.1 ; ival1(11:end)-0.1]
% UB=[ival1(1:10)+0.1 ; ival1(11:end)+0.1];

fixed_point_same_Z(ival1',fp,Fext,tF)

options = optimset('Display','iter','UseParallel',true','TolFun',1e-4,'TolX',1e-10,'MaxFunEvals',30000);
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(para) fixed_point_same_Z(para,fp,Fext,tF),ival1',[],[],[],[],LB,UB,@(para) compute_nonlcon_same_Z(para,fp,Fext,tF),options)

cost_verify=fixed_point_same_Z(x,fp,Fext,tF)
compute_nonlcon_same_Z(x,fp,Fext,tF)

x_save = [fp(1:5) ; x(1:5)' ; fp(11:end) ;x(6:end)'];

save('fixedpointforfivelink/Switching_Control/fixed_point_F=5_same_Z.mat','x_save');
