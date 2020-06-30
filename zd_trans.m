% D_bar=subs(D,[q6 q7 dq6 dq7],[c1-q1 c2 -dq1 0]);
% C_bar=subs(C,[q6 q7 dq6 dq7],[c1-q1 c2 -dq1 0])
D_bar=subs(D,[q6 q7 dq6 dq7],[c1 c2 0 0]);
C_bar=subs(C,[q6 q7 dq6 dq7],[c1 c2 0 0]);

Dl=D_bar(1:5,1:5);
Cl=C_bar(1:5,1:5);
Dml11=D_bar(1,6)
Cml11=C_bar(1,6)
ql=[q1;q2;q3;q4;q5];
dql=[dq1;dq2;dq3;dq4;dq5];
simplify(transpose(dql)*jacobian(Dl(1,:),ql)*dql-Cl(1,:)*dql)
check1=-dq1*jacobian(Dml11,ql)*dql;
check2=(1/2)*transpose(dql)*diff(Dl,q1)*dql;
dif=check1+check2
