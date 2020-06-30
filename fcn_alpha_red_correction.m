function [alpha,beta,theta_minus,theta_plus]=fcn_alpha_red_correction(ival,x_plus,N)
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

tau_plus = (c*x_plus(1:5) - theta_plus)/(theta_minus - theta_plus);
tau_dot_plus = (c*x_plus(6:10))/(theta_minus - theta_plus);

%==========================================================================
% Return to Z with a correction term
%==========================================================================
tau_return = 0.9;
hd = zeros(4,1);
for k=0:1:M
    hd = hd + alpha(:,k+1)*(factorial(M)/(factorial(k)*factorial(M-k)))*(tau_plus^k)*(1-tau_plus)^(M-k);
end
y_i = x_plus(2:5)-hd;

a= zeros(4,1);
for k=1:1:M-1
    a = a + (alpha(:,k+1))*(factorial(M)/(factorial(k)*factorial(M-k)))*(k*(tau_plus^(k-1))*(1-tau_plus)^(M-k) ...
        - (M-k)*(tau_plus^k)*(1-tau_plus)^(M-k-1));
end
y_dot_plus = x_plus(7:10) - (a - M*(alpha(:,1)*(1-tau_plus)^(M-1) - alpha(:,M+1)*(tau_plus^(M-1))))*tau_dot_plus;

clear a b
b = [y_i ; y_dot_plus/tau_dot_plus ; zeros(4,1) ; zeros(4,1) ; zeros(4,1) ; zeros(4,1)];

a = zeros(size(b,1),4*(N+1));
for k=0:1:N
    a(1:4,4*k+1:4*k+4) = (tau_plus^k)*eye(4);
    a(5:8,4*k+1:4*k+4) = k*(tau_plus^(k-1))*eye(4);
    a(9:12,4*k+1:4*k+4) = (tau_return^k)*eye(4);
    a(13:16,4*k+1:4*k+4) = k*(tau_return^(k-1))*eye(4);
    a(17:20,4*k+1:4*k+4) = k*(k-1)*(tau_return^(k-2))*eye(4);
    a(21:24,4*k+1:4*k+4) = k*(k-1)*(k-2)*(tau_return^(k-3))*eye(4);
end
beta_int = a\b;
beta = zeros(4,N+1);
for k=0:1:N
    beta(:,k+1) = beta_int(4*k+1:4*k+4);
end
% % % % for k=2:1:M-1
% % % %     a = a + (alpha(:,k+2)-alpha(:,k+1))*(factorial(M)*(tau_plus^k)*(1-tau_plus)^(M-1-k))...
% % % %         /(factorial(k)*factorial(M-1-k));
% % % % end
% % % % 
% % % % b_int=zeros(4,1);
% % % % for k=2:1:M
% % % %     b_int = b_int + alpha(:,k+1)*(factorial(M)*(tau_plus^k)*(1-tau_plus)^(M-k))...
% % % %         /(factorial(k)*factorial(M-k));
% % % % end
% % % % b = (M/(1-tau_plus))*(x_plus(2:5) - b_int);
% % % % 
% % % % alpha(:,2) = ( (x_plus(7:10)/tau_dot_plus) - alpha(:,3)*M*(M-1)*tau_plus*(1-tau_plus)^(M-2) - a + b)...
% % % %     /(M*((1-tau_plus)^(M-1)+tau_plus*(1-tau_plus)^(M-2)));

% for k=2:1:M-1
%     a = a + (alpha(:,k+1))*(factorial(M)/(factorial(k)*factorial(M-k)))*(k*(tau_plus^(k-1))*(1-tau_plus)^(M-k) ...
%         - (M-k)*(tau_plus^k)*(1-tau_plus)^(M-k-1));
% end
% 
% alpha(:,2) = ( (x_plus(7:10)/tau_dot_plus) + M*(alpha(:,1)*(1-tau_plus)^(M-1) - alpha(:,M+1)*(tau_plus^(M-1))) - a)...
%     /(M*((1-tau_plus)^(M-1) - (M-1)*tau_plus*(1-tau_plus)^(M-2)));
% 
% alpha_1_int = zeros(4,1);
% for k=1:1:M
%     alpha_1_int = alpha_1_int + alpha(:,k+1)*(factorial(M)*(tau_plus^k)*(1-tau_plus)^(M-k)...
%         /(factorial(k)*factorial(M-k)));
% end
% 
% alpha(:,1) = (x_plus(2:5) - alpha_1_int)/((1-tau_plus)^M);
% 
% hd = zeros(4,1);
% for k=0:1:M
%     hd = hd + alpha(:,k+1)*(factorial(M)/(factorial(k)*factorial(M-k)))*(tau_plus^k)*(1-tau_plus)^(M-k);
% end
% y_check = x_plus(2:5)-hd
% 
% for k=1:1:M-1
%     a = a + (alpha(:,k+1))*(factorial(M)/(factorial(k)*factorial(M-k)))*(k*(tau_plus^(k-1))*(1-tau_plus)^(M-k) ...
%         - (M-k)*(tau_plus^k)*(1-tau_plus)^(M-k-1));
% end
% a
% y_dot_check = x_plus(7:10) - (a - M*(alpha(:,1)*(1-tau_plus)^(M-1) - alpha(:,M+1)*(tau_plus^(M-1))))*tau_dot_plus
% tau_dot_plus

%==========================================================================
% Return to Z by the end of the step
%==========================================================================

% a = [-M*(1-tau_plus)^(M-1)*eye(4) (M*((1-tau_plus)^(M-1) - (M-1)*tau_plus*(1-tau_plus)^(M-2)))*eye(4) ;...
%      ((1-tau_plus)^M)*eye(4)       M*tau_plus*(1-tau_plus)^(M-1)*eye(4)];
%  
% b_1 = zeros(4,1);
% for k=2:1:M-1
%     b_1 = b_1 + (alpha(:,k+1))*(factorial(M)/(factorial(k)*factorial(M-k)))*(k*(tau_plus^(k-1))*(1-tau_plus)^(M-k) ...
%         - (M-k)*(tau_plus^k)*(1-tau_plus)^(M-k-1));
% end
% b = zeros(8,1);
% b(1:4,1) = ( (x_plus(7:10)/tau_dot_plus) - M*alpha(:,M+1)*(tau_plus^(M-1)) - b_1);
% 
% b_2 = zeros(4,1);
% for k=2:1:M
%     b_2 = b_2 + alpha(:,k+1)*(factorial(M)*(tau_plus^k)*(1-tau_plus)^(M-k)...
%         /(factorial(k)*factorial(M-k)));
% end
% b(5:8,1) = (x_plus(2:5) - b_2);
% 
% alpha_0_1 = a\b;
% alpha(:,1) = alpha_0_1(1:4,1);
% alpha(:,2) = alpha_0_1(5:8,1);

%==========================================================================
% Return to Z by the middle of the step
%==========================================================================
% b_1 = zeros(4,1);
% for k=4:1:M
%     b_1 = b_1 + alpha(:,k+1)*(factorial(M)*(tau_plus^k)*(1-tau_plus)^(M-k)...
%         /(factorial(k)*factorial(M-k)));
% end
% b_1 = x_plus(2:5) - b_1;
% 
% b_2 = zeros(4,1);
% b_2_1 = zeros(4,1);
% for k=4:1:M
%     b_2_1 = b_2_1 + alpha(:,k+1)*(factorial(M)*(0.5^k)*(1-0.5)^(M-k)...
%         /(factorial(k)*factorial(M-k)));
% end
% b_2_2 = zeros(4,1);
% for k=0:1:M
%     b_2_2 = b_2_2 + alpha(:,k+1)*(factorial(M)*(0.5^k)*(1-0.5)^(M-k)...
%         /(factorial(k)*factorial(M-k)));
% end
% b_2 = b_2_2 - b_2_1;
% 
% b_3 = zeros(4,1);
% b_3_1 = zeros(4,1);
% for k=4:1:M-1
%     b_3_1 = b_3_1 + (alpha(:,k+1))*(factorial(M)/(factorial(k)*factorial(M-k)))*(k*(tau_plus^(k-1))*(1-tau_plus)^(M-k) ...
%         - (M-k)*(tau_plus^k)*(1-tau_plus)^(M-k-1));
% end
% b_3 = ( (x_plus(7:10)/tau_dot_plus) - M*alpha(:,M+1)*(tau_plus^(M-1)) - b_3_1);
% 
% b_4 = zeros(4,1);
% b_4_1 = zeros(4,1);
% for k=4:1:M-1
%     b_4_1 = b_4_1 + (alpha(:,k+1))*(factorial(M)/(factorial(k)*factorial(M-k)))*(k*(0.5^(k-1))*(1-0.5)^(M-k) ...
%         - (M-k)*(0.5^k)*(1-0.5)^(M-k-1));
% end
% b_4_2 = zeros(4,1);
% for k=1:1:M-1
%     b_4_2 = b_4_2 + (alpha(:,k+1))*(factorial(M)/(factorial(k)*factorial(M-k)))*(k*(0.5^(k-1))*(1-0.5)^(M-k) ...
%         - (M-k)*(0.5^k)*(1-0.5)^(M-k-1));
% end
% b_4 = ( b_4_2 + M*(1-0.5)^(M-1)*alpha(:,1) - M*alpha(:,M+1)*(0.5^(M-1)) - b_4_1);
% 
% b = [b_1 ; b_2 ; b_3 ; b_4];
% 
% a=zeros(16);
% for k=0:1:3
%     a(1:4,4*k+1:4*k+4) = (factorial(M)*(tau_plus^k)*(1-tau_plus)^(M-k)/(factorial(k)*factorial(M-k)))*eye(4);
%     a(5:8,4*k+1:4*k+4) = (factorial(M)*(0.5^k)*(1-0.5)^(M-k)/(factorial(k)*factorial(M-k)))*eye(4);
%     if k==0
%         a(9:12,1:4) = -M*(1-tau_plus)^(M-1)*eye(4);
%         a(13:16,1:4) = -M*(1-0.5)^(M-1)*eye(4);
%     else
%         a(9:12,4*k+1:4*k+4) = (factorial(M)/(factorial(k)*factorial(M-k)))*(k*(tau_plus^(k-1))*(1-tau_plus)^(M-k) ...
%         - (M-k)*(tau_plus^k)*(1-tau_plus)^(M-k-1))*eye(4);
%     
%         a(13:16,4*k+1:4*k+4) = (factorial(M)/(factorial(k)*factorial(M-k)))*(k*(0.5^(k-1))*(1-0.5)^(M-k) ...
%         - (M-k)*(0.5^k)*(1-0.5)^(M-k-1))*eye(4);
%     end
% end
% alpha_0_1_2_3 = a\b;
% alpha(:,1) = alpha_0_1_2_3(1:4,1);
% alpha(:,2) = alpha_0_1_2_3(5:8,1);
% alpha(:,3) = alpha_0_1_2_3(9:12,1);
% alpha(:,4) = alpha_0_1_2_3(13:16,1);


%==========================================================================
% Display For Checking
%==========================================================================

% hd = zeros(4,1);
% for k=0:1:M
%     hd = hd + alpha(:,k+1)*(factorial(M)/(factorial(k)*factorial(M-k)))*(tau_plus^k)*(1-tau_plus)^(M-k);
% end
% y_check = x_plus(2:5)-hd
% 
% a= zeros(4,1);
% for k=1:1:M-1
%     a = a + (alpha(:,k+1))*(factorial(M)/(factorial(k)*factorial(M-k)))*(k*(tau_plus^(k-1))*(1-tau_plus)^(M-k) ...
%         - (M-k)*(tau_plus^k)*(1-tau_plus)^(M-k-1));
% end
% y_dot_check = x_plus(7:10) - (a - M*(alpha(:,1)*(1-tau_plus)^(M-1) - alpha(:,M+1)*(tau_plus^(M-1))))*tau_dot_plus
end