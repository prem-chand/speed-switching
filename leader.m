function [pL,pL_dot] = leader(t,x)
%LEADER Summary of this function goes here
%   Detailed explanation goes here
% pE = fcn_position_head(x(1:5));
global vL
% vL = 0.5;
pL = [vL*t-0.1131 ; 1.4066];
pL_dot = [vL; 0];

end

