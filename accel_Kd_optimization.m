clear;
clc;

M=6;
clear  theta_minus theta_plus betta XX event

LB=[-100 -100];
UB=[100 100];

% Genetic algorithm approach
mutationRate = 0.1;
options = gaoptimset('MutationFcn', {@mutationadaptfeasible, mutationRate},'TolFun',1e-2,'TolCon',1e-2,'generations',50,...
                                     'plotfcns',{@gaplotbestf,@gaplotbestindiv},'PopInitRange',[],'PopulationSize',10,...
                                      'FitnessLimit',1e-2,'InitialPopulation',[40 0]);
                                    
%
%[x,fval,exitflag,output,population,scores]=ga(@(para) fixed_point2(para,Fext,tF),26,[],[],[],[],LB,UB,@(para) compute_nonlcon(para,Fext,tF)...
 %                                             ,[],options)
[x,fval,exitflag,output,population,scores]=ga(@(para) KpKd_cost(para),2,[],[],[],[],LB,UB,[],options) %without constraints


% Fmincon approach
%[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(para) fixed_point2(para,Fext,tF),ival-0.1,[],[],[],[],LB,UB)
%[x,fval,exitflag,output,grad,hessian] = fminunc(@(para) fixed_point2(para,Fext,tF),ival)


%save('/Users/mohamadshafiee/Desktop/Reasearch/programming/5 link RABBIT/fixedpointforfivelink/fixed3_0N','x');

%save('/Users/mohamadshafiee/Desktop/Reasearch/programming/5 link RABBIT/fixedpointforfivelink/fixed_contrial','x');
