function [psi_H1,psi_H2,psi_V1,psi_V2] = fcn_GRF_QP_ineq(x,u,LgLfh_inv)
u=u;
dq = [x(6:10) ; 0 ; 0];
B = [zeros(1,4); eye(4) ; zeros(2,4) ];
[D,C,G,J,~] = fcn_extended_DCG(x);
psi_F1 = (J*(D\J'))\(J*(D\(C*dq+G-B*u)));
psi_F2 = -(J*(D\J'))\(J*(D\(B*LgLfh_inv)));
psi_H1 = psi_F1(1,:);
psi_H2 = psi_F2(1,:);
psi_V1 = psi_F1(2,:);
psi_V2 = psi_F2(2,:);
% GRF=GRF';
end

