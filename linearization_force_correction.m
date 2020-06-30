% function [A] = linearization_force_correction(N,F)
%LINEARIZATION_FORCE_CORRECTION Summary of this function goes here
%   Detailed explanation goes here
warning('OFF','all');
clear all
N=8;F=0 ;
fix=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0.mat');
% ival=fix.fpStates
ival=fix.x';
x_fix = ival(1:10);

clear  theta_minus theta_plus betta event XX M_bez
global theta_minus theta_plus betta event M_bez

[betta,beta_c,theta_minus,theta_plus] = fcn_alpha_red_correction(ival,x_fix,N);
beta_s = fcn_beta_correction(x_fix,x_fix,zeros(8,1));
%event=x_fix(2);
M_bez=6;
tF=[0 50];
Fext=[0  0];

% Integrator Constants
tfin = 5; relTol = 1e-5; absTol = 1e-4;
options   = odeset('Events',@(t,q) touchdown5_steps(t,q,0),...
    'RelTol',relTol,...
    'AbsTol',absTol,...
   'NormControl','off');
x = 1e-4; % perturbation

for i = 1:10
    x_fix = ival(1:10);
    %forward perturbation
    x_fix(i) = x_fix(i) + x;
    y0 = x_fix;
    [betta,beta_c,theta_minus,theta_plus] = fcn_alpha_red_correction(ival,x_fix,N);
    [~,~,~,ye,~]=...
        ode45(@(t,Z) return_map5(t,Z,Fext,tF,betta,beta_c,beta_s,theta_minus,theta_plus,M_bez,N),[0,tfin],y0,options);
    y0_plus = impact_map5(ye);
    
    %backward perturbation
    x_fix(i) = x_fix(i)-2*x;
    y0 = x_fix;
    [betta,beta_c,theta_minus,theta_plus] = fcn_alpha_red_correction(ival,x_fix,N);
    [~,~,~,ye,~]=...
        ode45(@(t,Z) return_map5(t,Z,Fext,tF,betta,beta_c,beta_s,theta_minus,theta_plus,M_bez,N),[0,tfin],y0,options);
    y0_minus = impact_map5(ye);
    
    % A is the Jacobian matrix of poincare map at fixed point, compute the
    % i-th column of A
    A(1:10,i) = (y0_plus-y0_minus)/(2*x);
    i
end
[~,v]=eig(A);
diag(v)
det(A-eye(10))


% % % % clear H
% % % % dF = 0.01;
% % % % for i = 1
% % % %     x_fix = ival(1:10);
% % % %     F_pert = F;
% % % %     %forward perturbation
% % % % %     x_fix(i) = x_fix(i) + x;
% % % %     F_pert = F(i) + dF;
% % % %     y0 = x_fix;
% % % %     beta_ctrl = zeros(8,1);
% % % %     beta_s_p = fcn_beta_correction(x_fix,x_fix,beta_ctrl);
% % % %     for k=1:1:4
% % % %         beta_s_p_tau(k,:) = fcn_coeffs_theta_to_tau(beta_s_p(k,:),theta_plus,theta_minus);
% % % %     end
% % % %     [betta,beta_c,theta_minus,theta_plus] = fcn_alpha_red_correction(ival,x_fix,N);
% % % %     [~,~,~,ye,~]=...
% % % %         ode45(@(t,Z) return_map5(t,Z,[F_pert;0],tF,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N),[0,tfin],y0,options);
% % % %     y0_plus = impact_map5(ye);
% % % %     
% % % %     %backward perturbation
% % % %     F_pert = F(i) - dF;
% % % % %     x_fix(i) = x_fix(i)-2*x;
% % % %     y0 = x_fix;
% % % %     beta_ctrl = zeros(8,1);
% % % %     beta_s_p = fcn_beta_correction(x_fix,x_fix,beta_ctrl);
% % % %     for k=1:1:4
% % % %         beta_s_p_tau(k,:) = fcn_coeffs_theta_to_tau(beta_s_p(k,:),theta_plus,theta_minus);
% % % %     end
% % % %     [betta,beta_c,theta_minus,theta_plus] = fcn_alpha_red_correction(ival,x_fix,N);
% % % %     [~,~,~,ye,~]=...
% % % %         ode45(@(t,Z) return_map5(t,Z,[F_pert;0],tF,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N),[0,tfin],y0,options);
% % % %     y0_minus = impact_map5(ye);
% % % %     
% % % %     % A is the Jacobian matrix of poincare map at fixed point, compute the
% % % %     % i-th column of A
% % % %     H(1:10,i) = (y0_plus-y0_minus)/(2*dF);
% % % % end
% % % % H
% % % % 
% % % % clear G beta_ctrl beta_s_p beta_s_p_tau
% % % % dbeta = 1e-2;
% % % % for i = 1:1:8
% % % %     x_fix = ival(1:10);
% % % %     %forward perturbation
% % % % %     x_fix(i) = x_fix(i) + x;
% % % %     y0 = x_fix;
% % % %     [betta,beta_c,theta_minus,theta_plus] = fcn_alpha_red_correction(ival,x_fix,N);
% % % %     beta_ctrl = zeros(8,1);
% % % %     beta_ctrl(i) = beta_ctrl(i) + dbeta;
% % % %     beta_s_p = fcn_beta_correction(x_fix,x_fix,beta_ctrl);
% % % %     for k=1:1:4
% % % %         beta_s_p_tau(k,:) = fcn_coeffs_theta_to_tau(beta_s_p(k,:),theta_plus,theta_minus);
% % % %     end
% % % %     beta_s_p_tau;
% % % %     [~,~,~,ye,~]=...
% % % %         ode45(@(t,Z) return_map5(t,Z,[F;0],tF,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N),[0,tfin],y0,options);
% % % %     y0_plus = impact_map5(ye);
% % % %     
% % % %     pause(1)
% % % %     %backward perturbation
% % % % %     x_fix(i) = x_fix(i)-2*x;
% % % %     y0 = x_fix;
% % % %     [betta,beta_c,theta_minus,theta_plus] = fcn_alpha_red_correction(ival,x_fix,N);
% % % %     beta_ctrl = zeros(12,1);
% % % %     beta_ctrl(i) = beta_ctrl(i) - dbeta;
% % % %     beta_s_p = fcn_beta_correction(x_fix,x_fix,beta_ctrl);
% % % %     for k=1:1:4
% % % %         beta_s_p_tau(k,:) = fcn_coeffs_theta_to_tau(beta_s_p(k,:),theta_plus,theta_minus);
% % % %     end
% % % %     [~,~,~,ye,~]=...
% % % %         ode45(@(t,Z) return_map5(t,Z,[F;0],tF,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N),[0,tfin],y0,options);
% % % %     y0_minus = impact_map5(ye);
% % % %     
% % % %     % A is the Jacobian matrix of poincare map at fixed point, compute the
% % % %     % i-th column of A
% % % %     G(1:10,i) = (y0_plus-y0_minus)/(2*dbeta);
% % % % end
% % % % G
% end

