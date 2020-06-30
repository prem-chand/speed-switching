% clear all;
clc;
syms q1 q2 q3 q4 q5 q6 q7 dq1 dq2 dq3 dq4 dq5 dq6 dq7 g lT lf lt MT Mf Mt IT If It pMT pMf pMt pMua pMfa lua Mua ...
    psh real Iua Mfa Ifa
q=[q1;q2;q3;q4;q5;q6;q7];
dq=[dq1; dq2; dq3;dq4;dq5;dq6;dq7];
g=9.8;lT=0.63;lf=0.4;lt=0.4;lua=0.25;lfa=0.25;
MT=12;Mf=6.8;Mt=3.2;Mua=1.36;Mfa=1;
IT=1.33;If=0.47;It=0.2;Iua=0.04;Ifa=0.03;
pMT=0.24;pMf=0.11;pMt=0.24;psh=0.38;pMua=0.125;pMfa=0.125;

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
KE=simple(KE)


%potential energy of links

PE=g*(Mf*(p_2(2)+p_3(2))+Mt*(p_4(2)+p_5(2))+MT*p_1(2)+Mua*p_6(2)+Mfa*p_7(2));
PE=simple(PE);

% gravity vector
Gl=transpose(jacobian(PE,q));
Gl=simple(Gl)

% mass-inertial matrix
%D=simple(jacobian(KE,dq).');
Dl=transpose(jacobian(KE,dq));
Dl=simple(jacobian(D,dq))

% Coriolis and centrigugal matrix

syms Cl real 
n=max(size(q));
for k=1:n
    for j=1:n
        Cl(k,j)=0*g;
        for i=1:n
            Cl(k,j)=Cl(k,j)+1/2*(diff(Dl(k,j),q(i))+...
                diff(Dl(k,i),q(j))-...
                diff(Dl(i,j),q(k)))*dq(i);
        end
    end
end
Cl=simple(Cl)

% input matrix
%B=[eye(4);zeros(1,4)];


%Jacobian of force
%P_Fext=[-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2);lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)];  % at hip
%P_Fext=[-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)+lT*sin(q1) lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)+lT*cos(q1)];  %at head
p_Fext=[-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)+psh*sin(q1)+lua*sin(pi-q1-q6)+lfa*sin(q1+q6-q7);...
         lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)+psh*cos(q1)-lua*cos(pi-q1-q6)+lfa*cos(q1+q6-q7)];         %at hand

p_shoulder=[-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)+psh*sin(q1);...
             lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)+psh*cos(q1)];         %at shoulder


     
% J_h=jacobian(p_Fext,q)
% J_h=simple(J_h)
% 
% J_sh=jacobian(p_shoulder,q)
% J_sh=simple(J_sh)

%computing Jdot which is a 2*7 matrix and need to be reconstructed for
%differentiation
% % Jdot_column=jacobian(J_h',q)*dq;
% % Jdot_column=simple(Jdot_column);
% % Jdot_h(1,:)=Jdot_column(1:7);
% % Jdot_h(2,:)=Jdot_column(8:14);
% % 
% % Jdot_column2=jacobian(J_sh',q)*dq;
% % Jdot_column2=simple(Jdot_column2);
% % Jdot_sh(1,:)=Jdot_column2(1:7);
% % Jdot_sh(2,:)=Jdot_column2(8:14);

%-------------------------------- Write the results
list_q_e  = {'q1','q(1)'; 'q2','q(2)'; 'q3','q(3)'; 'q4','q(4)';'q5','q(5)'; 'q6','q(6)';'q7','q(7)'};
             %'dq1','q(6)';'dq2','q(7)';'dq3','q(8)';'dq4','q(9)';'dq5','q(10)';'dq6','q(11)';'dq7','q(12)'};


p = mfilename('fullpath');
output_dir_m = [p '/'];

if not(exist(output_dir_m))
  mkdir(output_dir_m)
end

write_fcn_m([output_dir_m 'fcn_dynamics_matrices2.m'],{'q'},[list_q_e],{Gl,'Gl';Dl,'Dl';Cl,'Cl'});
%write_fcn_m([output_dir_m 'fcn_extended_C.m'],{'q','dq'},[list_q_e;list_dq_e],{C_e,'C_e'});



