function [alpha,theta_minus,theta_plus]=fcn_alpha_red(ival)
XX=ival(1:10);
c=[1 1 0 0.5 0];

theta_minus=XX(1)+XX(3)+0.5*XX(5); %c*q-
theta_plus=c*XX(1:5);             %c*q+
thetadot_plus=c*XX(6:10);        %c*q.+
M=(max(size(ival))-10)/4+2;

alpha=zeros(4,M+1);

alpha(:,1)=XX(2:5);               %y=h0-hd   hd=b(s=0)=beta0  y=0---> beta0=h0=[q2 q3 q4 q5]
H0=[zeros(4,1) eye(4)];
alpha(:,2)=H0*XX(6:10)*(theta_minus-theta_plus)/thetadot_plus/M   +    alpha(:,1); % equation 6.28 book and 77 paper
alpha(:,end)=[XX(3);XX(2);XX(5);XX(4)];

alpha(1,3:M)=ival(11:8+M);
alpha(2,3:M)=ival(9+M:6+2*M);
alpha(3,3:M)=ival(7+2*M:4+3*M);
alpha(4,3:M)=ival(5+3*M:2+4*M);

end