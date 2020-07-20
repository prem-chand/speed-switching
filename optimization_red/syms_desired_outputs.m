%function generate_stance_controller

clear all
clc

syms Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y9 Y10 theta real

syms s1 s2 s3 s4 theta_plus theta_minus
syms a01 a11 a21 a31 a41 a51 a61 
syms a03 a13 a23 a33 a43 a53 a63
syms a04 a14 a24 a34 a44 a54 a64
syms a02 a12 a22 a32 a42 a52 a62 

% -- Bezier coefficients
a1 = [a01 ; a11 ; a21 ; a31 ; a41 ; a51 ; a61];
a2 = [a02 ; a12 ; a22 ; a32 ; a42 ; a52 ; a62];
a3 = [a03 ; a13 ; a23 ; a33 ; a43 ; a53 ; a63];
a4 = [a04 ; a14 ; a24 ; a34 ; a44 ; a54 ; a64];


% -- Stance states
q_s  = [Y1 Y2 Y3 Y4 Y5];
dq_s = [Y6 Y7 Y8 Y9 Y10];

% -- Controlled variabes
controlled_vars = [Y2;Y3;Y4;Y5];

% -- Stance phase monotonic variable
% theta=[1 1 0 0.5 0]*q_s';

%==========================================================================
% -- Stance Controller 2: (qmLS->Bezier,qTor->Bezier)
%==========================================================================
% This controller is valid for theta_1<=theta<=theta_2

%--------------------------------------------------------------------------
% -- First output Bezier wrt s1
%--------------------------------------------------------------------------
[hd1, dhd1]=poly_eval_bezier(s1,6,a1);

% Dhd1_Ds1 = jacobian(hd1,s1);
% Dhd1_Ds1 = simple(Dhd1_Ds1);

% s1       =(theta-theta_plus)/(theta_minus-theta_plus);
% Ds1_Dq   = jacobian(s1,q_s); %the reason we dont take derivative wrt dq_s is that it will return 5 zeros which will later be multiplied by...
%                              %the last 5 rows of f(x).
% 
% Dhd1_Dq  = Dhd1_Ds1 * Ds1_Dq;  %chain rule
% Dhd1_Dq  = simple(Dhd1_Dq);
% Dh1_Dq   = jacobian(controlled_vars(1),q_s') - Dhd1_Dq;
% Dh1_Dq   = simple(Dh1_Dq);


%--------------------------------------------------------------------------
% -- Second output Bezier-part wrt s2
%--------------------------------------------------------------------------

[hd2,dhd2] = poly_eval_bezier(s2,6,a2);

% Dhd2_Ds2 = jacobian(hd2,s2);
% Dhd2_Ds2 = simple(Dhd2_Ds2);

% s2       = (theta - theta_plus) / (theta_minus - theta_plus);
% Ds2_Dq   = jacobian(s2,q_s);
% 
% Dhd2_Dq  = Dhd2_Ds2 * Ds2_Dq;  %chain rule
% Dhd2_Dq  = simple(Dhd2_Dq);
% Dh2_Dq   = jacobian(controlled_vars(2),q_s') - Dhd2_Dq;
% Dh2_Dq   = simple(Dh2_Dq);
%--------------------------------------------------------------------------
% -- Third output Bezier-part wrt s3
%--------------------------------------------------------------------------

[hd3,dhd3] = poly_eval_bezier(s3,6,a3);

% Dhd3_Ds3 = jacobian(hd3,s3);
% Dhd3_Ds3 = simple(Dhd3_Ds3);

% s3      = (theta - theta_plus) / (theta_minus - theta_plus);
% Ds3_Dq   = jacobian(s3,q_s);
% 
% Dhd3_Dq  = Dhd3_Ds3 * Ds3_Dq;  %chain rule
% Dhd3_Dq  = simple(Dhd3_Dq);
% Dh3_Dq   = jacobian(controlled_vars(3),q_s') - Dhd3_Dq;
% Dh3_Dq   = simple(Dh3_Dq);

%--------------------------------------------------------------------------
% -- Fourth output Bezier-part wrt s4
%--------------------------------------------------------------------------

[hd4,dhd4] = poly_eval_bezier(s4,6,a4);

% Dhd4_Ds4 = jacobian(hd4,s4);
% Dhd4_Ds4 = simple(Dhd4_Ds4);

% s4      = (theta - theta_plus) / (theta_minus - theta_plus);
% Ds4_Dq   = jacobian(s4,q_s);
% 
% Dhd4_Dq  = Dhd4_Ds4 * Ds4_Dq;  %chain rule
% Dhd4_Dq  = simple(Dhd4_Dq);
% Dh4_Dq   = jacobian(controlled_vars(4),q_s') - Dhd4_Dq;
% Dh4_Dq   = simple(Dh4_Dq);

hd = [subs(hd1, s1, (theta - theta_plus) / (theta_minus - theta_plus)),...
      subs(hd2, s2, (theta - theta_plus) / (theta_minus - theta_plus)),...
      subs(hd3, s3, (theta - theta_plus) / (theta_minus - theta_plus)),...
      subs(hd4, s4, (theta - theta_plus) / (theta_minus - theta_plus))];
hd = simplify(hd);

dhd = [subs(dhd1, s1, (theta - theta_plus) / (theta_minus - theta_plus)),...
      subs(dhd2, s2, (theta - theta_plus) / (theta_minus - theta_plus)),...
      subs(dhd3, s3, (theta - theta_plus) / (theta_minus - theta_plus)),...
      subs(dhd4, s4, (theta - theta_plus) / (theta_minus - theta_plus))]/(theta_minus - theta_plus);
dhd = simplify(dhd);

matlabFunction(hd,dhd,'file','optimization_red/desired_outputs.m','Optimize',false,...
    'Vars',{[a01, a11, a21, a31, a41, a51, a61;...
    a02, a12, a22, a32, a42, a52, a62;...
    a03, a13, a23, a33, a43, a53, a63;...
    a04, a14, a24, a34, a44, a54, a64],...
    theta, theta_minus, theta_plus})



































%--------------------------------------------------------------------------
% -- Combine both outputs
%--------------------------------------------------------------------------

% -- First order Lie derivative (combine both outputs)
% Dh_Dq  = [Dh1_Dq ; Dh2_Dq ; Dh3_Dq ; Dh4_Dq];
% Dh_Dq  = simple(Dh_Dq);
% 
% Lfh    = Dh_Dq * dq_s';
% Lfh    = simple(Lfh);
% hs=[hd1;hd2;hd3;hd4]
% h=[subs(Y2-hd1,{'s1'},{'(Y4/2 + Y1 + Y2 - theta_plus) / (theta_minus - theta_plus)'});
%    subs(Y3-hd2,{'s2'},{'(Y4/2 + Y1 + Y2 - theta_plus) / (theta_minus - theta_plus)'});
%    subs(Y4-hd3,{'s3'},{'(Y4/2 + Y1 + Y2 - theta_plus) / (theta_minus - theta_plus)'});
%    subs(Y5-hd4,{'s4'},{'(Y4/2 + Y1 + Y2 - theta_plus) / (theta_minus - theta_plus)'})]
% % -- Derivative of the first order Lie derivative wrt the states 
% Lfh_sub = subs(Lfh,{'s1','s2','s3','s4'},...
%     {'(Y4/2 + Y1 + Y2 - theta_plus) / (theta_minus - theta_plus)',...
%      '(Y4/2 + Y1 + Y2 - theta_plus) / (theta_minus - theta_plus)',...
%      '(Y4/2 + Y1 + Y2 - theta_plus) / (theta_minus - theta_plus)',...
%      '(Y4/2 + Y1 + Y2 - theta_plus) / (theta_minus - theta_plus)'})
%  
% Lfh=Lfh_sub;
% DLfh_Dx = jacobian(Lfh_sub,[q_s dq_s]');
% DLfh_Dx=simple(DLfh_Dx)
% 
% %-------------------------------- Write the results
% list_q_e  = {'(theta_plus - theta_minus)','deltheta';'Y1','Y(1)'; 'Y2','Y(2)'; 'Y3','Y(3)'; 'Y4','Y(4)'; ...
%     'Y5','Y(5)'; 'Y6','Y(6)';'Y7','Y(7)';'Y8','Y(8)';'Y9','Y(9)';'Y10','Y(10)';'Y(1)0','Y(10)';...
%     'a01','alpha(1,1)';'a11','alpha(1,2)';'a21','alpha(1,3)';'a31','alpha(1,4)';'a41','alpha(1,5)';'a51','alpha(1,6)';'a61','alpha(1,7)';...
%     'a02','alpha(2,1)';'a12','alpha(2,2)';'a22','alpha(2,3)';'a32','alpha(2,4)';'a42','alpha(2,5)';'a52','alpha(2,6)';'a62','alpha(2,7)';...
%     'a03','alpha(3,1)';'a13','alpha(3,2)';'a23','alpha(3,3)';'a33','alpha(3,4)';'a43','alpha(3,5)';'a53','alpha(3,6)';'a63','alpha(3,7)';...
%     'a04','alpha(4,1)';'a14','alpha(4,2)';'a24','alpha(4,3)';'a34','alpha(4,4)';'a44','alpha(4,5)';'a54','alpha(4,6)';'a64','alpha(4,7)';};
% %list_dq_e  = {'dqLA','dq(1)'; 'dqLS','dq(2)'; 'dqmLS','dq(3)'; ...
% %    'dqTor','dq(4)'; 'dxCOM','dq(5)'; 'dyCOM','dq(6)';};
% 
% p = mfilename('fullpath');
% output_dir_m = [p '/'];
% 
% if not(exist(output_dir_m))
%   mkdir(output_dir_m)
% end
% 
% write_fcn_m([output_dir_m 'fcn_liederivative.m'],{'Y','alpha'},[list_q_e],{h,'h';Lfh,'Lfh';DLfh_Dx,'DLfh_Dx';Dh_Dq,'Dh_Dq'});
