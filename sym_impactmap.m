clear
syms q1 q2 q3 q4 q5 dq1 dq2 dq3 dq4 dq5 real %g lT lf lt MT Mf Mt IT If It pMT pMf pMt real
syms z1 z2 dz1 dz2 real
g=9.8;lT=0.63;lf=0.4;lt=0.4;MT=12;Mf=6.8;Mt=3.2;IT=1.33;If=0.47;It=0.2;pMT=0.24;pMf=0.11;pMt=0.24;
q=[q1;q2;q3;q4;q5;z1;z2];
dq=[dq1; dq2; dq3;dq4;dq5;dz1;dz2];

p_toe=[z1;z2];

p_4=[z1-(lt-pMt)*sin(pi-q1-q2-q4);z2+(lt-pMt)*cos(pi-q1-q2-q4)];
v_4=jacobian(p_4,q)*dq;

p_2=[z1-lt*sin(pi-q1-q2-q4)-(lf-pMf)*sin(pi-q1-q2);z2+lt*cos(pi-q1-q2-q4)+(lf-pMf)*cos(pi-q1-q2)];
v_2=jacobian(p_2,q)*dq;

p_1=[z1-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)+pMT*sin(q1);z2+lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)+pMT*cos(q1)];
v_1=jacobian(p_1,q)*dq;

p_3=[z1-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)-pMf*sin(q1+q3-pi);z2+lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)-pMf*cos(q1+q3-pi)];
v_3=jacobian(p_3,q)*dq;

p_5=[z1-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)-lf*sin(q1+q3-pi)-pMt*sin(q1+q3+q5-pi);z2+lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)-lf*cos(q1+q3-pi)-pMt*cos(q1+q3+q5-pi)];
v_5=jacobian(p_5,q)*dq;


% kinetic energy of links
KE_1=1/2*MT*v_1.'*v_1+0.5*IT*dq1^2;
KE_2=1/2*Mf*v_2.'*v_2+0.5*If*(dq1+dq2)^2;
KE_3=1/2*Mf*v_3.'*v_3+0.5*If*(dq1+dq3)^2;
KE_4=1/2*Mt*v_4.'*v_4+0.5*It*(dq1+dq2+dq4)^2;
KE_5=1/2*Mt*v_5.'*v_5+0.5*It*(dq1+dq3+dq5)^2;

% total kinetic energy
KE=KE_1+KE_2+KE_3+KE_4+KE_5;
KE=simple(KE);

%potential energy of links
PE=g*(Mf*(p_2(2)+p_3(2))+Mt*(p_4(2)+p_5(2))+MT*p_1(2));
PE=simple(PE);

De=simple(jacobian(KE,dq).');
De=simple(jacobian(De,dq))


% gravity vector
Ge=jacobian(PE,q)';
Ge=simple(Ge)


% Coriolis and centrigugal matrix

syms C real
n=max(size(q));
for k=1:n
    for j=1:n
        C(k,j)=0;
        for i=1:n
            C(k,j)=C(k,j)+1/2*(diff(De(k,j),q(i))+...
                diff(De(k,i),q(j))-...
                diff(De(i,j),q(k)))*dq(i);
        end
    end
end
Ce=simple(C)

% foot jacobian
jacF  = jacobian(p_toe,q);
% Derivative of foot jacobian
reshape_jacF = reshape(jacF,2*7,1)
for i = 1:length(reshape_jacF);
    djacF(i) = jacobian(reshape_jacF(i),q)*dq;
end
djacF = reshape(djacF,2,7)


gama=[z1-lt*sin(pi-q1-q2-q4)-lf*sin(pi-q1-q2)-lf*sin(q1+q3-pi)-lt*sin(q1+q3+q5-pi);z2+lt*cos(pi-q1-q2-q4)+lf*cos(pi-q1-q2)-lf*cos(q1+q3-pi)-lt*cos(q1+q3+q5-pi)];
E=simple(jacobian(gama,q))

%-------------------------------- Write the results
list_q_e  = {'q1','q(1)'; 'q2','q(2)'; 'q3','q(3)'; 'q4','q(4)';'q5','q(5)'};


p = mfilename('fullpath');
output_dir_m = [p '/'];

if not(exist(output_dir_m))
  mkdir(output_dir_m)
end

%write_fcn_m([output_dir_m 'fcn_extended_D_E.m'],{'q'},[list_q_e],{De,'De';E,'E'});
write_fcn_m([output_dir_m 'fcn_extended_DCG.m'],{'q'},[list_q_e],{De,'De';Ce,'Ce';Ge,'Ge';jacF,'jacF';djacF,'djacF'});

%write_fcn_m([output_dir_m 'fcn_extended_C.m'],{'q','dq'},[list_q_e;list_dq_e],{C_e,'C_e'});


