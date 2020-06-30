clear all;
syms c l theta alpha theta1 theta2 l1 l2 ymx ymy dymx dymy q1 dq1 dl dtheta beta
q = [q1 ; ymx ; ymy];
dq = [dq1 ; dymx ; dymy];
% % q = [q1 ; l ; theta];
% % dq = [dq1 ; dl ; dtheta];
% % % expression=(c*sin(alpha)+l*sin(theta+alpha)-l1*sin(theta1))^2+...
% % %            (c*cos(alpha)+l*cos(theta+alpha)-l1*cos(theta1))^2;
% % % expression=simplify(expression);
% % A=c+l*cos(theta);
% % B=l*sin(theta);
% % delta=acos(A/sqrt(A^2+B^2));
% % % theta1=alpha+delta-acos((l1^2-l2^2+c^2+l^2+2*c*l*cos(theta))/(2*l1*sqrt(A^2+B^2)));
% % q6=beta+delta-acos((l1^2-l2^2+c^2+l^2+2*c*l*cos(theta))/(2*l1*sqrt(A^2+B^2)));
% % % theta2=theta1-atan((c*sin(alpha)+l*sin(theta+alpha)-l1*sin(theta1))/(c*cos(alpha)+l*cos(theta+alpha)-l1*cos(theta1)));
% % q7=q6-atan((c*sin(beta)+l*sin(theta+beta)-l1*sin(q6))/(c*cos(beta)+l*cos(theta+beta)-l1*cos(q6)));
% % dq6=jacobian(q6,q)*dq
% % dq7=jacobian(q7,q)*dq

% theta1=simplify(subs(expand(theta1),[l theta],[sqrt(yx^2+yy^2) atan(yx/yy)-alpha]))
% theta2=simplify(subs(expand(theta2),[l theta],[sqrt(yx^2+yy^2) atan(yx/yy)-alpha]))
% % x=[l ; theta];
% % dx=[dl ; dtheta];
% % qm=[theta1 ; theta2];
% % dqm=jacobian(qm,x)*dx;
% % simplify(dqm)


delta = acos((c+ymx*sin(q1+beta)+ymy*cos(q1+beta))/(sqrt(c^2+ymx^2+ymy^2+2*c*(ymx*sin(q1+beta)+ymy*cos(q1+beta)))));
q6 = beta + delta - acos((sqrt(c^2+ymx^2+ymy^2+2*c*(ymx*sin(q1+beta)+ymy*cos(q1+beta))))/(2*l1));
q6 = simplify(q6)
q7 = q1+q6-atan((c*sin(q1+beta)+ymx-l1*sin(q1+q6))/(c*cos(q1+beta)+ymy-l1*cos(q1+q6)));
q7 = simplify(q7)
dq6 = jacobian(q6,q)*dq;
dq6 = simplify(dq6)
dq7 = jacobian(q7,q)*dq;
dq7 = simplify(dq7)

