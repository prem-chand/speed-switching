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
num_steps = 1;
start = 79;
path = start;
warning('off','all');

clear  theta_minus theta_plus betta XX event
global theta_minus theta_plus betta event M_bez

fix=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0.mat');
ival=fix.x;
XX=ival(1:10)';

% % % % edge = load('./Fixed_Point_Lib/fp_F=0/edge_coarse_CLF_108_wts.mat');
% % % % edge_matrix_coarse = edge.edge_matrix;
% % % % gr = digraph(edge_matrix_coarse);
% % % % bins = conncomp(gr);
% % % % 
% % % % fps = load('./Fixed_Point_Lib/fp_F=0/fps_lib_coarse.mat');
% % % % data = fps.data;
% % % % XX = data(start,3:12)';
% % % % desired_speed = [0.42 0.81];
% % % % for i = 1:1:length(desired_speed)
% % % %     [~,desired_point(i)] = min(abs(data(:,1)-desired_speed(i)));
% % % %     [path_intermediate,cost(i)] = shortestpath(gr,path(end),desired_point(i));
% % % %     path = [path path_intermediate(2:end)];
% % % % end
% % % % size(path)
% % % % k = start;
% % % % zeta = data(start,13);

% % % jacob=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0_time_beta_jacobian.mat');
% % % G = jacob.G;

jacob=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0_time_beta_jacobian_m8.mat');
G = jacob.G;


% path = shortestpath(graph,79,39);
% path = [start 28];
% num_steps = 20*length(path);

% index = [0.43 0.7 0.75 0.8];

M_bez=6;
% XX=[];
fp_original = ival(1:10)';
c=[1 1 0 0.5 0];
dwell_time = 0;
last_step = 0;
count = 2;
delta_z = 0.9159;
epsilon = 2;
i = 1;
count2 = 1;
% for i=1:num_steps
t_history = [];
u_history = [];
% while count~=-1
k=1;
for i=1:1:num_steps
    clear s theta_minus theta_plus betta event

% % % %       if i >= last_step + dwell_time
% % % %           if count>length(path)
% % % %               break;
% % % %           elseif path(count) == desired_point(count2)
% % % %               epsilon = 0.5;
% % % %               count2 = count2 + 1;
% % % %           else
% % % %               epsilon = 2;
% % % %           end
% % % %           epsilon
% % % %           dwell_time = ceil(0.5*(log((abs(zeta-data(path(count),13))/epsilon)+1))/(log(1/delta_z)))
% % % %           k = path(count);
% % % %           last_step = i;
% % % %           count = count+1;
% % % %       end
% % % %           
%     k = path(i);
% % % %     desired_speed = data(k,1);
%     [betta,beta_c,theta_minus,theta_plus]=fcn_alpha_red_correction(ival',XX(:,i),N);
    [betta,beta_c,theta_minus,theta_plus]=fcn_alpha_red_correction(ival',fp_original,N);
%     beta_ctrl = pinv(G)*data(k,2);
    beta_ctrl = pinv(G)*-2.5e-2;
%     beta_ctrl = pinv(G)*0;
    beta_s_p = fcn_beta_correction_new(XX(:,1),XX(:,i),beta_ctrl);
    for j=1:1:4
        beta_s_p_tau(j,:) = fcn_coeffs_theta_to_tau_m8(beta_s_p(j,:),theta_plus,theta_minus);
    end
    Fext = [0 ; 0];
    options = odeset('Events',@(t,q)touchdown5_steps(t,q,0),'RelTol',1e-5,'AbsTol',1e-4);
    [t,y,te,ye,ie] = ode45(@(t,Z) return_map5_PD(t+tplot,Z,Fext,tF,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N),[0 20],XX(:,i),options);
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
    
    zeta = 0.5*(gamma0(ye)*ye(6:10)')^2;
    
    XX_pi = ye';
    [XX(:,i+1) , Fimp] = impact_map5(ye);
     gam0_plus=gamma0(XX(:,i+1));
    gam0_minus=gamma0(ye);
    delta_zero=(gam0_plus*XX(6:10,i+1))/(gam0_minus*ye(6:10)');
    timpact(i)=tplot;
    impactline=[timpact(i) timpact(i)];
    clear u h hdot hc hcdot h_original hdot_original GRF s friction_cone;
    clear output_error_d output_error_c output_error_s output_error
    for tp=1:size(t)
        [fx,gx,gforce]=fcn_stance_dynamics5(y(tp,:));
        u(tp,:)=fcn_stance_controller5_correction_PD(y(tp,:)',fx,gx,gforce,Fext,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N);
        if M_bez==6
            [h(tp,:),hdot(tp,:),~]=fcn_liederivative(y(tp,:),betta,theta_minus,theta_plus,M_bez);
            [hc(tp,:),hcdot(tp,:),~]=fcn_liederivative_correction(y(tp,:)',beta_c,theta_minus,theta_plus,N,0.5);
            [hs(tp,:),hsdot(tp,:),~]=fcn_liederivative_correction(y(tp,:)',beta_s_p_tau,theta_minus,theta_plus,8,0.9);
            [h_original(tp,:),hdot_original(tp,:)]=fcn_liederivative(y(tp,:),betta,theta_minus,theta_plus,M_bez);
        elseif M_bez==5
            [h(tp,:),hdot(tp,:)]=fcn_liederivative_m5(y(tp,:),betta);
        end
        s(tp)=(c*y(tp,1:5)'-theta_plus)/(theta_minus-theta_plus);
        GRF(tp,:)=fcn_GRF(y(tp,:),u(tp,:));
        friction_cone(tp) = abs(GRF(tp,1)/GRF(tp,2));
        cm(tp,:)=fcn_cm_kinematics(y(tp,1:5));
        theta(tp) = c*y(tp,1:5)';
    end
   
    
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
    average_speed(i)=step_length(i,1)/te;
    
    if i>=2
        accel(i)=(average_speed(i)-average_speed(i-1))/te;
    end
    
    %     friction_cone = [friction_cone abs(Fimp(1)/Fimp(2))];
    fric_max = max(friction_cone);
    GRF_hor_min = min(GRF(:,1));
    %     figure(1)
    %     plot(t+tplot,GRF(:,1));hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    [u_max,u_max_index] = max(max(abs(u')));
        u_max;
        desired_speed = 0;
    fprintf('i=%1.0f,\t desired_speed=%1.4f,\t average_speed=%1.4f,\t umax= %3.2f,\t friction_cone= %2.4f \t max_output_error= %1.5f \n',i,desired_speed,average_speed(i),u_max,fric_max,max(max(output_error(:,:))));
    c*y(u_max_index,1:5)';
    t(u_max_index);
    %     error=norm(XX(:,i+1)-ival(1:10)')
    error=norm(XX(:,i+1)-XX(:,i));
    
    
%     figure(1)
%     subplot(2,2,1);plot(t+tplot,u(:,1));title('u1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%     subplot(2,2,2);plot(t+tplot,u(:,2));title('u2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%     subplot(2,2,3);plot(t+tplot,u(:,3));title('u3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%     subplot(2,2,4);plot(t+tplot,u(:,4));title('u4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    
    tplot=tplot+te;
    
        if i==num_steps
%             figure(1)
%             subplot(2,4,1);plot(t+tplot,output_error_c(:,1));title('h1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,2);plot(t+tplot,output_error_c(:,2));title('h2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,3);plot(t+tplot,output_error_c(:,3));title('h3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,4);plot(t+tplot,output_error_c(:,4));title('h4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%     
%             subplot(2,4,5);plot(t+tplot,output_error_c(:,5));title('h^.1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,6);plot(t+tplot,output_error_c(:,6));title('h^.2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,7);plot(t+tplot,output_error_c(:,7));title('h^.3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,8);plot(t+tplot,output_error_c(:,8));title('h^.4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    
            figure(2)
            subplot(2,4,1);plot(t+tplot,output_error_s(:,1));title('h1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
            subplot(2,4,2);plot(t+tplot,output_error_s(:,2));title('h2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
            subplot(2,4,3);plot(t+tplot,output_error_s(:,3));title('h3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
            subplot(2,4,4);plot(t+tplot,output_error_s(:,4));title('h4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    
            subplot(2,4,5);plot(t+tplot,output_error_s(:,5));title('h^.1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
            subplot(2,4,6);plot(t+tplot,output_error_s(:,6));title('h^.2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
            subplot(2,4,7);plot(t+tplot,output_error_s(:,7));title('h^.3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
            subplot(2,4,8);plot(t+tplot,output_error_s(:,8));title('h^.4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%     
%             figure(3)
%             subplot(2,4,1);plot(t+tplot,output_error_d(:,1));title('h1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,2);plot(t+tplot,output_error_d(:,2));title('h2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,3);plot(t+tplot,output_error_d(:,3));title('h3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,4);plot(t+tplot,output_error_d(:,4));title('h4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%     
%             subplot(2,4,5);plot(t+tplot,output_error_d(:,5));title('h^.1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,6);plot(t+tplot,output_error_d(:,6));title('h^.2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,7);plot(t+tplot,output_error_d(:,7));title('h^.3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,8);plot(t+tplot,output_error_d(:,8));title('h^.4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%     
%             figure(4)
%             subplot(2,4,1);plot(t+tplot,output_error(:,1));title('h1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,2);plot(t+tplot,output_error(:,2));title('h2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,3);plot(t+tplot,output_error(:,3));title('h3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,4);plot(t+tplot,output_error(:,4));title('h4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%     
%             subplot(2,4,5);plot(t+tplot,output_error(:,5));title('h^.1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,6);plot(t+tplot,output_error(:,6));title('h^.2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,7);plot(t+tplot,output_error(:,7));title('h^.3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%             subplot(2,4,8);plot(t+tplot,output_error(:,8));title('h^.4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        end
% break;
        hs_save = hs;
        theta_save = theta;
        t_history = [t_history ; t+tplot];
        u_history = [u_history ; u];
    clear output_error norm_output_error hs hsdot theta;
    i = i+1;
end
% % figure(5)
% % for i = 1:1:size(edge_matrix_coarse,1)
% %     for j = 1:1:size(edge_matrix_coarse,2)
% %         if edge_matrix_coarse(i,j)~=0
% %             plot(data(j,1),data(i,1),'b*','MarkerSize',1);
% %             hold on;
% %         end
% %     end
% % end
% % plot([data(1,1) data(end,1)],[data(1,1) data(end,1)],'r');
% % hold off;
% XX;
% timpact(i+1)=tplot;
