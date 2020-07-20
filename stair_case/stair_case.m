clc; clear all; close all

num_steps = 10;
stair_ht = randi([-2,3],num_steps,1)/100;
ctrl_seq = randi([-2,3],num_steps,1)/100;

file = 'fixedpointforfivelink/Switching_Control/fixed_point_d=0cm.mat';
fix = open(file);
ival = fix.x;
x0 = ival(1:10)';

XX = [];
XX(:,1) = x0;
for i = 1:num_steps
    XX(:,i+1) = Poincare_stair(XX(:,i),stair_ht(i),ctrl_seq(i));
    
    
    
end