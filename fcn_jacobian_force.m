function J = fcn_jacobian_force(q)
%FCN_JACOBIAN_FORCE Summary of this function goes here
%   Detailed explanation goes here

J(1,1)=(63*cos(q(1)))/100 - (2*cos(q(1) + q(2)))/5 - (2*cos(q(1) + q(2) + q(4)))/5;
  J(1,2)=- (2*cos(q(1) + q(2) + q(4)))/5 - (2*cos(q(1) + q(2)))/5;
  J(1,3)=0;
  J(1,4)=-(2*cos(q(1) + q(2) + q(4)))/5;
  J(1,5)=0;
  J(2,1)=(2*sin(q(1) + q(2) + q(4)))/5 + (2*sin(q(1) + q(2)))/5 - (63*sin(q(1)))/100;
  J(2,2)=(2*sin(q(1) + q(2) + q(4)))/5 + (2*sin(q(1) + q(2)))/5;
  J(2,3)=0;
  J(2,4)=(2*sin(q(1) + q(2) + q(4)))/5;
  J(2,5)=0;
end

