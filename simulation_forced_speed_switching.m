clear all;
close all;
clc;
M=6;
N = 8;
te=0;
tplot=0;timpact=0;
ntouch=0;ptouch=0;
accel(1)=0;
simout.t=[];simout.stk=[];simout.h=[]; simout.swk=[]; simout.swf=[]; simout.head=[];
tF=[0 500];
Fext=[];
angle = 0;
num_steps = 50;
start = 1;
vL = 0.82;

clear  theta_minus theta_plus betta XX event
% global theta_minus theta_plus betta event M_bez

fix=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0.mat');
ival=fix.x;
% XX=ival(1:10)'

% fps = load('./Fixed_Point_Lib/fp_F=0/fps_lib_coarse.mat');
fps = load('./Fixed_Point_Lib/fp_F=0/fps_lib_fine.mat');
data = fps.data;
XX = data(start,3:12)';

% edge = load('./Fixed_Point_Lib/fp_F=0/edge_coarse_CLF_108.mat');
% edge_matrix_coarse = edge.edge_matrix;
% graph = digraph(edge_matrix_coarse);
% bins = conncomp(graph);
%
% fps = load('./Fixed_Point_Lib/fp_F=0/fps_lib_coarse.mat');
% data = fps.data;
% XX = data(start,3:12)';
% desired_speed = [0.8 0.6 0.5];
% for i = 1:1:length(desired_speed)
%     [~,desired_point(i)] = min(abs(data(:,1)-desired_speed(i)));
%     path_intermediate = shortestpath(graph,path(end),desired_point(i));
%     path = [path path_intermediate(2:end)];
% end
% path
% k = start;
% zeta = data(start,13);

jacob=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0_time_beta_jacobian.mat');
G = jacob.G;

M_bez=6;
% XX=[];
% k=1;
c=[1 1 0 0.5 0];
% speed_ctrl = [7.4e-3 6.5e-3];
k=start;
stance_foot_zero = 0;
for i=1:num_steps
    clear s theta_minus theta_plus betta event
    
    desired_speed = data(k,1);
%     [betta,beta_c,theta_minus,theta_plus]=fcn_alpha_red_correction(ival',XX(:,1),N);
        [betta,beta_c,theta_minus,theta_plus]=fcn_alpha_red_correction(ival',XX(:,i),N);
    beta_ctrl = pinv(G)*data(k,2);
    beta_s_p = fcn_beta_correction(XX(:,1),XX(:,i),beta_ctrl);
    for j=1:1:4
        beta_s_p_tau(j,:) = fcn_coeffs_theta_to_tau(beta_s_p(j,:),theta_plus,theta_minus);
    end
    Fext = [0 ; 0];
    options = odeset('Events',@(t,q)touchdown5_steps(t,q,0),'RelTol',1e-5,'AbsTol',1e-4);
    [t,y,te,ye,ie] = ode45(@(t,Z) return_map5_cooperation(t+tplot,Z,stance_foot_zero,tF,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N),[0 20],XX(:,i),options);
    simout.t=[simout.t;t+tplot];
    [pstk ph pswk pswf phead vh Menergy]=out_kinematics(y);
    
    simout.stk = [simout.stk  ; pstk];
    simout.h   = [simout.h    ; ph];
    simout.swk = [simout.swk  ; pswk];
    simout.swf = [simout.swf  ; pswf];
    simout.head= [simout.head ; phead];
    
    ntouch=ntouch+size(pswf,1);
    ptouch=ptouch+pswf(size(pswf,1),1);
    str_lng(i,1)=ntouch;
    str_lng(i,2)=ptouch;
    
    XX_pi = ye';
    [XX(:,i+1) , Fimp] = impact_map5(ye);
    gam0_plus=gamma0(XX(:,i+1));
    gam0_minus=gamma0(ye);
    delta_zero=(gam0_plus*XX(6:10,i+1))/(gam0_minus*ye(6:10)');
    timpact(i)=tplot;
    impactline=[timpact(i) timpact(i)];
    
    clear u h hdot hc hcdot h_original hdot_original GRF s;
    clear k1 k2 k3 Fvir zeta_f Vzero Fext_plot pL pE Fext_norm friction_cone theta
    for tp=1:size(t)
        [fx,gx,gforce]=fcn_stance_dynamics5(y(tp,:));
        beta_ctrl = pinv(G)*data(k,2);
        beta_s_p = fcn_beta_correction(XX(:,1),XX(:,i),beta_ctrl);
        for j=1:1:4
            beta_s_p_tau(j,:) = fcn_coeffs_theta_to_tau(beta_s_p(j,:),theta_plus,theta_minus);
        end
        
        u(tp,:)=fcn_stance_controller5_correction(y(tp,:)',fx,gx,gforce,Fext,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N);
        if M_bez==6
            [h(tp,:),hdot(tp,:),~]=fcn_liederivative(y(tp,:),betta,theta_minus,theta_plus,M_bez);
            [hc(tp,:),hcdot(tp,:),~]=fcn_liederivative_correction(y(tp,:)',beta_c,theta_minus,theta_plus,N,0.5);
            [hs(tp,:),hsdot(tp,:),~]=fcn_liederivative_correction(y(tp,:)',beta_s_p_tau,theta_minus,theta_plus,7,0.9);
            [h_original(tp,:),hdot_original(tp,:)]=fcn_liederivative(y(tp,:),betta,theta_minus,theta_plus,M_bez);
        elseif M_bez==5
            [h(tp,:),hdot(tp,:)]=fcn_liederivative_m5(y(tp,:),betta);
        end
        s(tp)=(c*y(tp,1:5)'-theta_plus)/(theta_minus-theta_plus);
        GRF(tp,:)=fcn_GRF(y(tp,:),u(tp,:));
        friction_cone(tp) = abs(GRF(tp,1)/GRF(tp,2));
        [pL(:,tp),pL_dot(:,tp)] = leader(t(tp)+tplot,y(tp,:));
        pE(tp,:) = fcn_position_head(y(tp,1:5)) + [stance_foot_zero 0];
        Fext_plot(:,tp) = compute_coop_force(y(tp,:)',pL(:,tp),pL_dot(:,tp),stance_foot_zero);
        Fext_norm(tp) = norm(Fext_plot(:,tp));
        diff(tp) = pL(1,tp)-pE(tp,1);
        
        cm(tp,:)=fcn_cm_kinematics(y(tp,1:5));
        theta(tp)=c*y(tp,1:5)';
        s(tp)=(theta(tp)-theta_plus)/(theta_minus-theta_plus);
        k2(tp)=kappa2(y(tp,1:5));
        k3(tp,:)=kappa3(y(tp,1:5));
        for j=1:1:size(data,1)
            beta_ctrl = pinv(G)*data(j,2);
            beta_s_p = fcn_beta_correction(XX(:,1),XX(:,i),beta_ctrl);
            beta_s_p_tau = fcn_hs_theta_to_tau(beta_s_p,theta_plus,theta_minus);
%             for m=1:1:4
%                 beta_s_p_tau2(m,:) = fcn_coeffs_theta_to_tau(beta_s_p(m,:),theta_plus,theta_minus);
%             end
            k1(j,tp)=c*([fcn_Dh_Dq(s(tp),betta,(theta_plus-theta_minus))-fcn_Dhc_Dq(y(tp,:)',beta_s_p_tau,theta_minus,theta_plus,7,0.9)-fcn_Dhc_Dq(y(tp,:)',beta_c,theta_minus,theta_plus,N,0.5) ; gamma0(y(tp,1:5))]\[zeros(4,1) ; 1]);
            Fvir(j,tp)=(k3(tp,:)*Fext_plot(:,tp))/k1(j,tp);
        end
        ksi1(tp)=theta(tp);
        ksi2(tp)=gamma0(y(tp,1:5)')*y(tp,6:10)';
    end
    
    ksi2_minus = gamma0(y(end,1:5)')*y(end,6:10)';
    zeta_minus = 0.5*ksi2_minus^2;
    
            for j=1:1:size(data,1)
                Wf(j)=trapz(theta,Fvir(j,:));
                zeta_f(j) = data(k,13) + Wf(j)/(1-delta_zero^2);
            end
    %     max(Wf)
    %     min(Wf)
    %     Wf
            [~,k] = min(abs(zeta_f'-data(:,13)));
    %         abs(zeta_f'-data(:,13))
    
    %     plot(t,diff);
    %     plot(t+tplot,pL(1,:));
    %     hold on
    %     plot(t+tplot,pE(:,1));
    % %    plot(t,Fext_plot(1,:));
    %     pause(0.1);
    
    output_error_c=[hc hcdot];
    output_error_s=[hs hsdot];
    output_error_d=[h hdot];
    output_error=[(h-hc-hs) (hdot-hcdot-hsdot)];
    norm_output_error_c = max(sqrt(sum(output_error_c.^2,2)));
    norm_output_error_s = max(sqrt(sum(output_error_s.^2,2)));
    norm_output_error_d = max(sqrt(sum(output_error_d.^2,2)));
    %     norm_output_error = max(sqrt(sum(output_error.^2,2)))
    norm_output_error = max(max(abs(output_error)));
    %     Fimp = [cosd(angle) -sind(angle) ; sind(angle) cosd(angle)]*Fimp';
    
    step_length(i,:)=fcn_position_swingfoot(ye);
    stance_foot_zero = step_length(i,1) + stance_foot_zero;
    average_speed(i)=step_length(i,1)/te;
    end_time(i)=te;
    
    if i>=2
        accel(i)=(average_speed(i)-average_speed(i-1))/te;
    end
    
    Force_RMS(i) = rms(Fext_norm);
    
    %     [~,k] = min(abs(data(:,1)-average_speed(i)));
%     [~,k] = min(abs(data(:,13)-zeta_minus));
    
    %     if i>1
    %         speed_check = average_speed(i) + accel(i)*te;
    %         [~,k] = min(abs(data(:,1)-speed_check));
    %     end
    
    
    
    %     friction_cone = [friction_cone abs(Fimp(1)/Fimp(2))];
    fric_max = max(friction_cone);
    GRF_hor_min = min(GRF(:,1));
    %     figure(1)
    %     plot(t+tplot,GRF(:,1));hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    [u_max,u_max_index] = max(max(abs(u')));
    u_max;
    c*y(u_max_index,1:5)';
    t(u_max_index);
    
    %     error=norm(XX(:,i+1)-ival(1:10)')
    error=norm(XX(:,i+1)-XX(:,i));
    fprintf('i=%1.0f,\t desired_speed=%1.4f,\t average_speed=%1.4f,\t umax= %3.2f,\t friction_cone= %2.4f \t max_output_error= %1.5f \n',i,desired_speed,average_speed(i),u_max,fric_max,max(max(output_error(:,:))));
    
    
    
    %     figure(1)
    %     subplot(2,2,1);plot(t+tplot,u(:,1));title('u1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %     subplot(2,2,2);plot(t+tplot,u(:,2));title('u2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %     subplot(2,2,3);plot(t+tplot,u(:,3));title('u3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %     subplot(2,2,4);plot(t+tplot,u(:,4));title('u4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    
    tplot=tplot+te;
    
    %     if i==num_steps
    %         figure(1)
    %         subplot(2,4,1);plot(t+tplot,output_error_c(:,1));title('h1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,2);plot(t+tplot,output_error_c(:,2));title('h2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,3);plot(t+tplot,output_error_c(:,3));title('h3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,4);plot(t+tplot,output_error_c(:,4));title('h4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %
    %         subplot(2,4,5);plot(t+tplot,output_error_c(:,5));title('h^.1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,6);plot(t+tplot,output_error_c(:,6));title('h^.2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,7);plot(t+tplot,output_error_c(:,7));title('h^.3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,8);plot(t+tplot,output_error_c(:,8));title('h^.4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %
    %         figure(2)
    %         subplot(2,4,1);plot(t+tplot,output_error_s(:,1));title('h1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,2);plot(t+tplot,output_error_s(:,2));title('h2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,3);plot(t+tplot,output_error_s(:,3));title('h3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,4);plot(t+tplot,output_error_s(:,4));title('h4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %
    %         subplot(2,4,5);plot(t+tplot,output_error_s(:,5));title('h^.1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,6);plot(t+tplot,output_error_s(:,6));title('h^.2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,7);plot(t+tplot,output_error_s(:,7));title('h^.3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,8);plot(t+tplot,output_error_s(:,8));title('h^.4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %
    %         figure(3)
    %         subplot(2,4,1);plot(t+tplot,output_error_d(:,1));title('h1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,2);plot(t+tplot,output_error_d(:,2));title('h2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,3);plot(t+tplot,output_error_d(:,3));title('h3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,4);plot(t+tplot,output_error_d(:,4));title('h4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %
    %         subplot(2,4,5);plot(t+tplot,output_error_d(:,5));title('h^.1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,6);plot(t+tplot,output_error_d(:,6));title('h^.2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,7);plot(t+tplot,output_error_d(:,7));title('h^.3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,8);plot(t+tplot,output_error_d(:,8));title('h^.4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %
    %         figure(4)
    %         subplot(2,4,1);plot(t+tplot,output_error(:,1));title('h1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,2);plot(t+tplot,output_error(:,2));title('h2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,3);plot(t+tplot,output_error(:,3));title('h3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,4);plot(t+tplot,output_error(:,4));title('h4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %
    %         subplot(2,4,5);plot(t+tplot,output_error(:,5));title('h^.1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,6);plot(t+tplot,output_error(:,6));title('h^.2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,7);plot(t+tplot,output_error(:,7));title('h^.3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %         subplot(2,4,8);plot(t+tplot,output_error(:,8));title('h^.4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %     end
    
    clear output_error norm_output_error hs hsdot;
end
plot(Force_RMS);
figure(2)
plot(average_speed)
hold on
plot(0:0.1:num_steps,vL*ones(size(0:0.1:num_steps,2)));

XX;
timpact(i+1)=tplot;
