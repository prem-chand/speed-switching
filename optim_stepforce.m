clear;
close all;
clc;
M=6;
te=0;
tplot=0;timpact=0;
ntouch=0;ptouch=0;
accel(1)=0;
simout.t=[];simout.stk=[];simout.h=[]; simout.swk=[]; simout.swf=[]; simout.head=[];

tF=[0 .2241 .2242 0.4379 0.4380  50];
Fext=[0 0    10   10     0      0];
    
fix=open('/Users/mohamadshafiee/Desktop/Reasearch/programming/5 link RABBIT/fixedpointforfivelink/fixed_point_9.mat');
ival=fix.fpStates';
fix=open('/Users/mohamadshafiee/Desktop/Reasearch/programming/5 link RABBIT/fixedpointforfivelink/converg_2.mat');
ival2=fix.x;
fix=open('/Users/mohamadshafiee/Desktop/Reasearch/programming/5 link RABBIT/fixedpointforfivelink/converg_3.mat');
ival3=fix.x;

clear  theta_minus theta_plus betta XX event XX_noforce
global theta_minus theta_plus betta event XX_noforce
XX=[];
XX=ival(1:10)';
event=XX(2);

c=[1 1 0 0.5 0];
[betta,theta_minus,theta_plus]=fcn_alpha_red(ival');


for i=1:3
clear s
theta_minus=XX(1,i)+XX(3,i)+0.5*XX(5,i); %c*q-
theta_plus=c*XX(1:5,i);                  %c*q+

%thetadot_plus=c*XX(6:10,i)        %c*q.+
%med=H0*XX(6:10,i)/thetadot_plus/M*(theta_minus-theta_plus)
%betta(:,2)=med(1:4)+betta(:,1)   %equation 6.16a

options = odeset('Events',@touchdown5,'RelTol',1e-5,'AbsTol',1e-5);
[t,y,te,ye,ie] = ode45(@(t,Z) return_map5(t+tplot,Z,Fext,tF),[0 20],XX(:,i),options);

XX(:,i+1)=impact_map5(ye);
tplot=tplot+te;
timpact(i)=tplot

end

XX_noforce(:,1)=XX(:,4);
XX_noforce(:,2)=XX(:,1);
XX;
%animate

%clear  theta_minus theta_plus betta XX event


%tF=[0  1 1.0001 8 8.001 20];
%Fext=[0 0 50    50 0      0];
%tF=  [0  0.2 0.2001 0.5 0.5001 20];
%Fext=[0   0    10    10   10     10];
tF=[0 10];
Fext=[0 0];
%LB=[];
%UB=[];
%ival(6:10)=ival(6:10)+0.1;
%LB=ival'-0.2;
%UB=ival'+0.2;

LB=ival(11:26)'-0.3;
UB=ival(11:26)'+0.3;


% Genetic algorithm approach
mutationRate = 0.2;
options = gaoptimset('MutationFcn', {@mutationadaptfeasible, mutationRate},'TolFun',1e-6,'TolCon',1e-6,'generations',400,...
                                     'plotfcns',{@gaplotbestf,@gaplotbestindiv},'PopInitRange',[],'PopulationSize',16,...
                                      'FitnessLimit',8e-4,'InitialPopulation',[ival(11:26);ival2;ival3]);
                               
                                    
%
%[x,fval,exitflag,output,population,scores]=ga(@(para) fixed_convergance_noforce(para,Fext,tF),16,[],[],[],[],LB,UB,@(para) compute_nonlcon(para,Fext,tF)...
%                                              ,[],options)
%[x,fval,exitflag,output,population,scores]=ga(@(para) fixed_convergance_noforce(para,Fext,tF),16,[],[],[],[],LB,UB,[],options) %without constraints

%save('/Users/mohamadshafiee/Desktop/Reasearch/programming/5 link RABBIT/fixedpointforfivelink/converg_4','x');

