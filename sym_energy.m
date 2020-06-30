clear
syms q1 q2 q3 q4 q5 dq1 dq2 dq3 dq4 dq5 g lT lf lt MT Mf Mt IT If It pMT pMf pMt real
q=[q1;q2;q3;q4;q5];
dq=[dq1; dq2; dq3;dq4;dq5];
g=9.8;lT=0.63;lf=0.4;lt=0.4;MT=12;Mf=6.8;Mt=3.2;IT=1.33;If=0.47;It=0.2;pMT=0.24;pMf=0.11;pMt=0.24;

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

% kinetic energy of links
KE_1=1/2*MT*v_1.'*v_1+0.5*IT*dq1^2;
KE_2=1/2*Mf*v_2.'*v_2+0.5*If*(dq1+dq2)^2;
KE_3=1/2*Mf*v_3.'*v_3+0.5*If*(dq1+dq3)^2;
KE_4=1/2*Mt*v_4.'*v_4+0.5*It*(dq1+dq2+dq4)^2;
KE_5=1/2*Mt*v_5.'*v_5+0.5*It*(dq1+dq3+dq5)^2;

% total kinetic energy
KE=KE_1+KE_2+KE_3+KE_4+KE_5;
KE=simple(KE)


%potential energy of links

PE=g*(Mf*(p_2(2)+p_3(2))+Mt*(p_4(2)+p_5(2))+MT*p_1(2));
PE=simple(PE);

ME=PE+KE;


%-------------------------------- Write the results
list_q_e  = {'dq1','q(6)';'dq2','q(7)';'dq3','q(8)';'dq4','q(9)';'dq5','q(10)';...
    'q1','q(1)'; 'q2','q(2)'; 'q3','q(3)'; 'q4','q(4)';'q5','q(5)'; };
%list_dq_e  = {'dqLA','dq(1)'; 'dqLS','dq(2)'; 'dqmLS','dq(3)'; ...
%    'dqTor','dq(4)'; 'dxCOM','dq(5)'; 'dyCOM','dq(6)';};

p = mfilename('fullpath');
output_dir_m = [p '/'];

if not(exist(output_dir_m))
  mkdir(output_dir_m)
end

write_fcn_m([output_dir_m 'fcn_energy.m'],{'q'},[list_q_e],{ME,'ME';KE,'KE';PE,'PE';});
%write_fcn_m([output_dir_m 'fcn_extended_C.m'],{'q','dq'},[list_q_e;list_dq_e],{C_e,'C_e'});



