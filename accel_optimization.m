clear;
clc;

M=6;
clear  theta_minus theta_plus betta XX event

LB=[0 0 0 0 0];
UB=[15 15 15 15 15];
%LB(6:10)=0.01;
%UB(6:10)=1.5;


LB
UB
% Genetic algorithm approach
mutationRate = 0.4;
options = gaoptimset('MutationFcn', {@mutationadaptfeasible, mutationRate},'TolFun',1e-2,'TolCon',1e-2,'generations',40,...
                                     'plotfcns',{@gaplotbestf,@gaplotbestindiv},'PopInitRange',[],'PopulationSize',10,...
                                      'FitnessLimit',1e-2,'InitialPopulation',[]);
                                    
%
%[x,fval,exitflag,output,population,scores]=ga(@(para) fixed_point2(para,Fext,tF),26,[],[],[],[],LB,UB,@(para) compute_nonlcon(para,Fext,tF)...
 %                                             ,[],options)
[x,fval,exitflag,output,population,scores]=ga(@(para) accel_cost(para),5,[],[],[],[],LB,UB,[],options) %without constraints


% Fmincon approach
%[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(para) fixed_point2(para,Fext,tF),ival-0.1,[],[],[],[],LB,UB)
%[x,fval,exitflag,output,grad,hessian] = fminunc(@(para) fixed_point2(para,Fext,tF),ival)


%save('/Users/mohamadshafiee/Desktop/Reasearch/programming/5 link RABBIT/fixedpointforfivelink/fixed3_0N','x');

%save('/Users/mohamadshafiee/Desktop/Reasearch/programming/5 link RABBIT/fixedpointforfivelink/fixed_contrial','x');
