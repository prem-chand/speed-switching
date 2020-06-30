clear all;
% clc;
syms q1 q2 q3 q4 q5 q6 q7 dq1 dq2 dq3 dq4 dq5 dq6 dq7 g lT lf lt MT Mf Mt IT If It pMT pMf pMt pMua pMfa lfa lua Mua ...
    psh real Iua Mfa Ifa c1 c2
q=[q1;q2;q3;q4;q5;q6;q7];
dq=[dq1; dq2; dq3;dq4;dq5;dq6;dq7];
% g=9.8;lT=0.63;lf=0.4;lt=0.4;lua=0.25;lfa=0.25;
% MT=12;Mf=6.8;Mt=3.2;Mua=1.36;Mfa=1;
% IT=1.33;If=0.47;It=0.2;Iua=0.04;Ifa=0.03;
% pMT=0.24;pMf=0.11;pMt=0.24;psh=0.38;pMua=0.125;pMfa=0.125;

%position of center of mass of each link as a function of states
p_4=[-(lt-pMt)*sin(pi-q1-q2-q4);(lt-pMt)*cos(pi-q1-q2-q4)];
v_4=jacobian(p_4,q)*dq;

p_2=[-lt*sin(pi-q1-q2-q4)-(lf-pMf)*sin(pi-q1-q2);lt*cos(pi-q1-q2-q4)+(lf-pMf)*cos(pi-q1-q2)];
v_2=jacobian(p_2,q)*dq;

p_1=[-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)+pMT*sin(q1);lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)+pMT*cos(q1)];
v_1=jacobian(p_1,q)*dq;

p_3=[-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)-pMf*sin(q1+q3-pi);lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)-pMf*cos(q1+q3-pi)];
v_3=jacobian(p_3,q)*dq;

p_5=[-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)-lf*sin(q1+q3-pi)-pMt*sin(q1+q3+q5-pi);lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)-lf*cos(q1+q3-pi)-pMt*cos(q1+q3+q5-pi)];
v_5=jacobian(p_5,q)*dq;

p_6=[-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)+psh*sin(q1)+pMua*sin(pi-q1-q6);lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)+psh*cos(q1)-pMua*cos(pi-q1-q6)];
v_6=jacobian(p_6,q)*dq;

p_7=[-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)+psh*sin(q1)+lua*sin(pi-q1-q6)+pMfa*sin(q1+q6-q7);lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)+psh*cos(q1)-lua*cos(pi-q1-q6)+pMfa*cos(q1+q6-q7)];
v_7=jacobian(p_7,q)*dq;

% kinetic energy of links
KE_1=1/2*MT*transpose(v_1)*v_1+0.5*IT*dq1^2;
KE_2=1/2*Mf*transpose(v_2)*v_2+0.5*If*(dq1+dq2)^2;
KE_3=1/2*Mf*transpose(v_3)*v_3+0.5*If*(dq1+dq3)^2;
KE_4=1/2*Mt*transpose(v_4)*v_4+0.5*It*(dq1+dq2+dq4)^2;
KE_5=1/2*Mt*transpose(v_5)*v_5+0.5*It*(dq1+dq3+dq5)^2;

KE_6=1/2*Mua*transpose(v_6)*v_6+0.5*Iua*(dq1+dq6)^2;
KE_7=1/2*Mfa*transpose(v_7)*v_7+0.5*Ifa*(dq1+dq6-dq7)^2;

% total kinetic energy
KE=KE_1+KE_2+KE_3+KE_4+KE_5+KE_6+KE_7;
KE=simplify(KE)


%potential energy of links

PE=g*(Mf*(p_2(2)+p_3(2))+Mt*(p_4(2)+p_5(2))+MT*p_1(2)+Mua*p_6(2)+Mfa*p_7(2));
PE=simple(PE);

% gravity vector
G=jacobian(PE,q)';
G=simplify(G)

% mass-inertial matrix
%D=simple(jacobian(KE,dq).');
D=jacobian(KE,dq).';
D=simplify(jacobian(D,dq))

% Coriolis and centrigugal matrix

syms C real
n=max(size(q));
for k=1:n
    for j=1:n
        C(k,j)=0*g;
        for i=1:n
            C(k,j)=C(k,j)+1/2*(diff(D(k,j),q(i))+...
                diff(D(k,i),q(j))-...
                diff(D(i,j),q(k)))*dq(i);
        end
    end
end
C=simplify(C)


% D_bar=subs(D,[q6 q7 dq6 dq7],[c1-q1 c2 -dq1 0]);
% C_bar=subs(C,[q6 q7 dq6 dq7],[c1-q1 c2 -dq1 0])
% Dl=D_bar(1:5,1:5);
% Dml11=D_bar(1,6);
% ql=[q1;q2;q3;q4;q5];
% dql=[dq1;dq2;dq3;dq4;dq5];
% check1=-dq1*jacobian(Dml11,ql)*dql
% check2=(1/2)*transpose(dql)*diff(Dl,q1)*dql
% dif=check1+check2
% Cml11=C_bar(1,6)



% Dl=[D(1:5,1)-D(1:5,6) D(1:5,2:5)];
% Dl=simplify(Dl)
% input matrix
%B=[eye(4);zeros(1,4)];


%Jacobian of force
%P_Fext=[-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2);lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)];  % at hip
%P_Fext=[-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)+lT*sin(q1) lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)+lT*cos(q1)];  %at head
p_Fext=[-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)+psh*sin(q1)+lua*sin(pi-q1-q6)+lfa*sin(q1+q6-q7);...
         lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)+psh*cos(q1)-lua*cos(pi-q1-q6)+lfa*cos(q1+q6-q7)];         %at hand

p_shoulder=[-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)+psh*sin(q1);...
             lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)+psh*cos(q1)];         %at shoulder


ym=p_Fext-p_shoulder;
Jm=jacobian(ym,q);
     
% J_h=jacobian(p_Fext,[q6 ; q7]);
% J_h=simple(J_h);
% 
% J_sh=jacobian(p_shoulder,[q6 ; q7]);
% J_sh=simple(J_sh);

% syms Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y9 Y10 real
syms dqLA dqLS dqmLS dqTor real

syms s1 s2 s3 s4 theta_plus theta_minus
syms a01  a11 a21 a31 a41 a51 a61 a71
syms a03  a13 a23 a33 a43 a53 a63
syms a04  a14 a24 a34 a44 a54 a64
syms a_12 a02 a12 a22 a32 a42 a52 a62 a72 %"_12" means "minus1,2":a12~=a_12 


% -- Bezier coefficients
a1 = [a01 ; a11 ; a21 ; a31 ; a41 ; a51 ; a61];
a2 = [a02 ; a12 ; a22 ; a32 ; a42 ; a52 ; a62];
a3 = [a03 ; a13 ; a23 ; a33 ; a43 ; a53 ; a63];
a4 = [a04 ; a14 ; a24 ; a34 ; a44 ; a54 ; a64];


% -- Stance states
q_s  = [q1;q2;q3;q4;q5];
dq_s = [dq1;dq2;dq3;dq4;dq5];

% -- Controlled variabes
controlled_vars = [q2;q3;q4;q5];

% -- Stance phase monotonic variable
theta=[1 1 0 0.5 0]*q_s;

%==========================================================================
% -- Stance Controller 2: (qmLS->Bezier,qTor->Bezier)
%==========================================================================
% This controller is valid for theta_1<=theta<=theta_2

%--------------------------------------------------------------------------
% -- First output Bezier wrt s1
%--------------------------------------------------------------------------
[hd1,dhd1]=poly_eval_bezier(((q1+q2+0.5*q4) - theta_plus) / (theta_minus - theta_plus),6,a1);

% Dhd1_Ds1 = jacobian(hd1,s1);
% Dhd1_Ds1 = simple(Dhd1_Ds1);
% 
% s1       =(theta-theta_plus)/(theta_minus-theta_plus);
% Ds1_Dq   = jacobian(s1,q_s); %the reason we dont take derivative wrt dq_s is that it will return 5 zeros which will later be multiplied by...
%                              %the last 5 rows of f(x).
% 
% Dhd1_Dq  = Dhd1_Ds1 * Ds1_Dq;  %chain rule
% Dhd1_Dq  = simple(Dhd1_Dq);
% Dh1_Dq   = jacobian(controlled_vars(1),q_s) - Dhd1_Dq;
% Dh1_Dq   = simple(Dh1_Dq);
h1=controlled_vars(1)-hd1;
h1=simplify(h1);
Dh1_Dq   = jacobian(h1,q_s);
Dh1_Dq   = simple(Dh1_Dq);


%--------------------------------------------------------------------------
% -- Second output Bezier-part wrt s2
%--------------------------------------------------------------------------

[hd2,dhd2] = poly_eval_bezier(((q1+q2+0.5*q4) - theta_plus) / (theta_minus - theta_plus),6,a2);

% Dhd2_Ds2 = jacobian(hd2,s2);
% Dhd2_Ds2 = simple(Dhd2_Ds2);
% 
% s2       = (theta - theta_plus) / (theta_minus - theta_plus);
% Ds2_Dq   = jacobian(s2,q_s);
% 
% Dhd2_Dq  = Dhd2_Ds2 * Ds2_Dq;  %chain rule
% Dhd2_Dq  = simple(Dhd2_Dq);
% Dh2_Dq   = jacobian(controlled_vars(2),q_s) - Dhd2_Dq;
% Dh2_Dq   = simple(Dh2_Dq);
h2=controlled_vars(2)-hd2;
h2=simplify(h2);
Dh2_Dq   = jacobian(h2,q_s);
Dh2_Dq   = simple(Dh2_Dq);
%--------------------------------------------------------------------------
% -- Third output Bezier-part wrt s3
%--------------------------------------------------------------------------

[hd3,dhd3] = poly_eval_bezier(((q1+q2+0.5*q4) - theta_plus) / (theta_minus - theta_plus),6,a3);

% Dhd3_Ds3 = jacobian(hd3,s3);
% Dhd3_Ds3 = simple(Dhd3_Ds3);
% 
% s3      = (theta - theta_plus) / (theta_minus - theta_plus);
% Ds3_Dq   = jacobian(s3,q_s);
% 
% Dhd3_Dq  = Dhd3_Ds3 * Ds3_Dq;  %chain rule
% Dhd3_Dq  = simple(Dhd3_Dq);
% Dh3_Dq   = jacobian(controlled_vars(3),q_s) - Dhd3_Dq;
% Dh3_Dq   = simple(Dh3_Dq);
h3=controlled_vars(3)-hd3;
h3=simplify(h3);
Dh3_Dq   = jacobian(h3,q_s);
Dh3_Dq   = simple(Dh3_Dq);

%--------------------------------------------------------------------------
% -- Fourth output Bezier-part wrt s4
%--------------------------------------------------------------------------

[hd4,dhd4] = poly_eval_bezier(((q1+q2+0.5*q4) - theta_plus) / (theta_minus - theta_plus),6,a4);

% Dhd4_Ds4 = jacobian(hd4,s4);
% Dhd4_Ds4 = simple(Dhd4_Ds4);
% 
% s4      = ((q1+q2+0.5*q4) - theta_plus) / (theta_minus - theta_plus);
% Ds4_Dq   = jacobian(s4,q_s);
% 
% Dhd4_Dq  = Dhd4_Ds4 * Ds4_Dq;  %chain rule
% Dhd4_Dq  = simple(Dhd4_Dq);
% Dh4_Dq   = jacobian(controlled_vars(4),q_s) - Dhd4_Dq;
% Dh4_Dq   = simple(Dh4_Dq);
h4=controlled_vars(4)-hd4;
h4=simplify(h4);
Dh4_Dq   = jacobian(h4,q_s);
Dh4_Dq   = simple(Dh4_Dq);
%--------------------------------------------------------------------------
% -- Combine both outputs
%--------------------------------------------------------------------------

% -- First order Lie derivative (combine both outputs)
Dh_Dq  = [Dh1_Dq ; Dh2_Dq ; Dh3_Dq ; Dh4_Dq];
Dh_Dq  = simple(Dh_Dq);

h=[h1 ; h2 ; h3 ; h4];

coord=[ h ; Dh_Dq*dq_s ; ym ; Jm*dq ; theta ; D(1,:)*dq]
coord_jacob=jacobian(coord,[q ; dq])
% determinant = det(coord_jacob)
% determinant = simplify(determinant)


