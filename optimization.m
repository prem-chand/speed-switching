clear;
clc;
warning('off','all');

M=6;
clear  theta_minus theta_plus betta XX
global theta_minus theta_plus betta XX

%tF=[0  1 1.0001 8 8.001 20];
%Fext=[0 0 50    50 0      0];
%tF=  [0  0.2 0.2001 0.5 0.5001 20];
%Fext=[0   0    10    10   10     10];
tF=[0 100];
Fext=[0 0];

fix=open('fixedpointforfivelink/fixed_point_9_gz_duanyi.mat');
ival=fix.x;
% ival=fix.fpStates';

% XX=[];
% XX=ival(1:10);

% c=[1 1 0 0.5 0];
% theta_minus=XX(1)+XX(3)+0.5*XX(5) %c*q-
% theta_plus=c*XX(1:5)             %c*q+
% thetadot_plus=c*XX(6:10);        %c*q.+
% betta=zeros(4,7);
% betta(:,1)=XX(2:5);               %y=h0-hd   hd=b(s=0)=beta0  y=0---> beta0=h0=[q2 q3 q4 q5]
% H=[zeros(4,1) eye(4);c];
% med=H*XX(6:10)/thetadot_plus/M*(theta_minus-theta_plus);
% betta(:,2)=med(1:4)+betta(:,1);   %equation 6.16a
% betta(:,7)=[XX(3);XX(2);XX(5);XX(4)];

LB=ival-0.5*abs(ival);
UB=ival+0.5*abs(ival);

options = optimset('Algorithm','sqp','display','iter','TolCon',1e-2,'MaxFunEvals',10000,'UseParallel',true,'TolFun',1e-5,'TolX',1e-8);

[x,fval,exitflag,output]  = fmincon(@(para) fixed_point2(para,Fext,tF),ival,[],[],[],[],LB,UB,[],options)

save('fixedpointforfivelink/fixed_point_9_sv','x');

%norm=fixed_point(ival(11:26),Fext,tF)
