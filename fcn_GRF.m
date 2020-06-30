function GRF = fcn_GRF(x,u)
u=u';
x=x';
dq = [x(6:10) ; 0 ; 0];
B = [zeros(1,4); eye(4) ; zeros(2,4) ];
[D,C,G,J,dJ] = fcn_extended_DCG(x);
GRF=(J*(D\J'))\(J*(D\(C*dq+G-B*u)));
% GRF=GRF';
end

