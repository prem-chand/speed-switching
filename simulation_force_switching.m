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
num_steps = 10;


% fix=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0.mat');
fix=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0_slow.mat');
ival=fix.x;
% XX=ival(1:10)'

jacob=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0_time_beta_jacobian.mat');
G = jacob.G;

pGH = 1.0e+06*[0.5783   -1.0944   -2.4011    4.5318   -0.6196    1.1696    1.7663   -3.3349]';

% hs_params = open('fixedpointforfivelink/Switching_Control/fixed_point_F=0_diff_speed2.mat');
% ival = hs_params.x_save(1:26)';
% beta_ctrl = hs_params.x_save(27:end)

clear  theta_minus theta_plus betta XX event
global theta_minus theta_plus betta event M_bez
M_bez=6;
XX=[];
XX=ival(1:10)';

c=[1 1 0 0.5 0];
% [betta,theta_minus,theta_plus]=fcn_alpha_red(ival');
% options = odeset('Events',@(t,q)touchdown5(t,q,angle),'RelTol',1e-5,'AbsTol',1e-4);
% [t,y,te,ye,ie] = ode45(@(t,Z) return_map5(t,Z,Fext,tF,betta,theta_minus,theta_plus,M_bez),[0 20],XX,options);
%
% zeta_fp = gamma0(y(end,:))*y(end,6:10)';
% d=0;
% F_actual = (ones(num_steps,1))*10;
% F_actual = (rand(num_steps,1))*50;
% F_actual = [0 0];
% F_actual = ones(num_steps,1)*-10;
% F = -10:5:10;
% XX_lib=[];
% for k=1:length(F)
%     fix_temp=open(strcat(strcat('fixedpointforfivelink/Switching_Control/fixed_point_F=',num2str(F(k))),'.mat'));
%     ival_temp = fix_temp.x;
%     XX_lib(:,k) = ival_temp(1:10)';
% end
for i=1:num_steps
    clear s theta_minus theta_plus betta event

    [betta,beta_c,theta_minus,theta_plus]=fcn_alpha_red_correction(ival',XX(:,i),N);
    beta_c = zeros([4 9]);
%     [betta,beta_c,theta_minus,theta_plus]=fcn_alpha_red_correction([XX(:,i)' ival(11:end)]',XX(:,i),N);
%     beta_ctrl = pinv(G)*(6e-3);
    beta_ctrl = -pGH*(0);
    beta_s_p = fcn_beta_correction(XX(:,1),XX(:,i),beta_ctrl);
    for k=1:1:4
        beta_s_p_tau(k,:) = fcn_coeffs_theta_to_tau(beta_s_p(k,:),theta_plus,theta_minus);
    end
    Fext = [5 ; 0];
%     Fext(:,i) = [F_actual(i) ; 0];
    options = odeset('Events',@(t,q)touchdown5_steps(t,q,0),'RelTol',1e-5,'AbsTol',1e-4);
    [t,y,te,ye,ie] = ode45(@(t,Z) return_map5(t+tplot,Z,Fext,tF,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N),[0 20],XX(:,i),options);
    simout.t=[simout.t;t+tplot];
%     te
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
    timpact(i)=tplot;
    impactline=[timpact(i) timpact(i)];
    
    clear u h hdot hc hcdot h_original hdot_original GRF s;
    for tp=1:size(t)
        [fx,gx,gforce]=fcn_stance_dynamics5(y(tp,:));
        %         Fext = force(t(tp)+tplot);
        %         Fext_plot = [interp1(tF,Fext,t(tp)+tplot);0];
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
        cm(tp,:)=fcn_cm_kinematics(y(tp,1:5));
    end
    %     section_state(:,i) = y(end,:)';
    %     V_eta(i) = [h(end,:),hdot(end,:)]*P*[h(end,:),hdot(end,:)]'
    %     zeta = gamma0(y(end,:))*y(end,6:10)';
    %     V_xi(i) = (zeta-zeta_fp)^2
    %     V(i) = 10*V_eta(i) + V_xi(i);
    
    output_error_c=[hc hcdot];
    output_error_s=[hs hsdot];
    output_error_d=[h hdot];
%     output_error=[h_original-(h-hc-hs) hdot_original-(hdot-hcdot-hsdot)];
    output_error=[(h-hc-hs) (hdot-hcdot-hsdot)];
    norm_output_error_c = max(sqrt(sum(output_error_c.^2,2)));
    norm_output_error_s = max(sqrt(sum(output_error_s.^2,2)));
    norm_output_error_d = max(sqrt(sum(output_error_d.^2,2)));
    norm_output_error = max(sqrt(sum(output_error.^2,2)));
    %     Fimp = [cosd(angle) -sind(angle) ; sind(angle) cosd(angle)]*Fimp';
    
    step_length(i,:)=fcn_position_swingfoot(ye);
    average_speed(i)=step_length(i,1)/te
    
    if i>=2
        accel(i)=(average_speed(i)-average_speed(i-1))/te;
    end
    
%     friction_cone = [friction_cone abs(Fimp(1)/Fimp(2))];
    fric_max = max(friction_cone)
    GRF_hor_min = min(GRF(:,1));
%     figure(1)
%     plot(t+tplot,GRF(:,1));hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    [u_max,u_max_index] = max(max(abs(u')));
    u_max
%     u_max
    c*y(u_max_index,1:5)';
    t(u_max_index);
    
    %     error=norm(XX(:,i+1)-ival(1:10)')
    error=norm(XX(:,i+1)-XX(:,i))
    
                
%          figure(2)
%         subplot(2,2,1);plot(t+tplot,u(:,1));title('u1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%         subplot(2,2,2);plot(t+tplot,u(:,2));title('u2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%         subplot(2,2,3);plot(t+tplot,u(:,3));title('u3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%         subplot(2,2,4);plot(t+tplot,u(:,4));title('u4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%     
    
        
%     figure(2)
%     subplot(2,4,1);plot(t+tplot,h_original(:,1));title('h1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%     subplot(2,4,2);plot(t+tplot,h_original(:,2));title('h2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%     subplot(2,4,3);plot(t+tplot,h_original(:,3));title('h3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%     subplot(2,4,4);plot(t+tplot,h_original(:,4));title('h4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%     
%     subplot(2,4,5);plot(t+tplot,hdot_original(:,1));title('h^.1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%     subplot(2,4,6);plot(t+tplot,hdot_original(:,2));title('h^.2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%     subplot(2,4,7);plot(t+tplot,hdot_original(:,3));title('h^.3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
%     subplot(2,4,8);plot(t+tplot,hdot_original(:,4));title('h^.4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    tplot=tplot+te;
    
    if i==num_steps
        figure(1)
        subplot(2,4,1);plot(t+tplot,output_error_c(:,1));title('h1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,2);plot(t+tplot,output_error_c(:,2));title('h2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,3);plot(t+tplot,output_error_c(:,3));title('h3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,4);plot(t+tplot,output_error_c(:,4));title('h4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    
        subplot(2,4,5);plot(t+tplot,output_error_c(:,5));title('h^.1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,6);plot(t+tplot,output_error_c(:,6));title('h^.2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,7);plot(t+tplot,output_error_c(:,7));title('h^.3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,8);plot(t+tplot,output_error_c(:,8));title('h^.4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        
        figure(2)
        subplot(2,4,1);plot(t+tplot,output_error_s(:,1));title('h1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,2);plot(t+tplot,output_error_s(:,2));title('h2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,3);plot(t+tplot,output_error_s(:,3));title('h3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,4);plot(t+tplot,output_error_s(:,4));title('h4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    
        subplot(2,4,5);plot(t+tplot,output_error_s(:,5));title('h^.1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,6);plot(t+tplot,output_error_s(:,6));title('h^.2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,7);plot(t+tplot,output_error_s(:,7));title('h^.3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,8);plot(t+tplot,output_error_s(:,8));title('h^.4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        
        figure(3)
        subplot(2,4,1);plot(t+tplot,output_error_d(:,1));title('h1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,2);plot(t+tplot,output_error_d(:,2));title('h2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,3);plot(t+tplot,output_error_d(:,3));title('h3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,4);plot(t+tplot,output_error_d(:,4));title('h4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    
        subplot(2,4,5);plot(t+tplot,output_error_d(:,5));title('h^.1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,6);plot(t+tplot,output_error_d(:,6));title('h^.2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,7);plot(t+tplot,output_error_d(:,7));title('h^.3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,8);plot(t+tplot,output_error_d(:,8));title('h^.4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        
        figure(4)
        subplot(2,4,1);plot(t+tplot,output_error(:,1));title('h1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,2);plot(t+tplot,output_error(:,2));title('h2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,3);plot(t+tplot,output_error(:,3));title('h3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,4);plot(t+tplot,output_error(:,4));title('h4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    
        subplot(2,4,5);plot(t+tplot,output_error(:,5));title('h^.1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,6);plot(t+tplot,output_error(:,6));title('h^.2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,7);plot(t+tplot,output_error(:,7));title('h^.3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        subplot(2,4,8);plot(t+tplot,output_error(:,8));title('h^.4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
        
        
        fp = XX(:,i+1);
%         save('fixedpointforfivelink/Switching_Control/fixed_point_F=0_0.8_beta.mat','betta');
%         save('fixedpointforfivelink/Switching_Control/fixed_point_F=0_0.8_beta_c.mat','beta_c');
%         save('fixedpointforfivelink/Switching_Control/fixed_point_F=0_0.8_beta_s.mat','beta_s_p_tau');
%         save('fixedpointforfivelink/Switching_Control/fixed_point_F=0_0.8_fp.mat','fp');
        beta_c
        beta_s_p_tau
        fp
    end
    clear output_error norm_output_error hs hsdot;
end
XX;
timpact(i+1)=tplot;
% figure(2)
% plot([1:num_steps],V_xi,'x','linewidth',2);
% hold on;
% plot([1:num_steps],V,'x','linewidth',2);
% xlabel('Step Number','fontsize',30);
% ylabel('V_x','fontsize',30);
% hold on;
% plot([1:num_steps],20*ones(1,num_steps),'--','color','k','linewidth',1);
% hold off;

% figure(2)
% plot([1:num_steps],V_eta,'x');
% figure(3)
% plot([1:num_steps],V_xi,'x');

% figure (9)
% %subplot(2,1,1);plot(timpact(2:end),average_speed,'rs');xlim([0 tplot]);xlabel('Time (s)','fontsize',16);ylabel('Average Speed (m/s)','fontsize',16);%title('Average Speed');
% %subplot(3,1,2);plot(timpact(2:end),accel,'--rs');xlim([0 tplot]);title('Average Accelaration');hold on;% stairs(timpact(2:end),a_des,'-ok');hold on
% %subplot(2,1,2);plot(tF,Fext,'--rs');ylim([0 12]);xlim([0 tplot]);xlabel('Time (s)','fontsize',16);ylabel('External Force','fontsize',16);%title('External Force');
% plot(timpact(2:end),average_speed,'rs');xlim([0 tplot]);xlabel('Time (s)','fontsize',16);ylabel('Average Speed (m/s)','fontsize',16);%title('Average Speed');
% figure (10)
% plot(tF,Fext,'--rs');ylim([0 12]);xlim([0 tplot]);xlabel('Time (s)','fontsize',16);ylabel('External Force','fontsize',16);%title('External Force');

mean_accel=(average_speed(i)-average_speed(1))/(tplot-timpact(2));
%sum(step_length(2:end,1));
%mean_accel=(sum(step_length(2:end,1))-average_speed(1)*(tplot-timpact(2)))*2/(tplot-timpact(2))^2

%figure(10)
%plot(timpact,XX(6,1:size(XX,2)-1));title('Angular Velocity of Torso after each impact');hold on;
%lim=ylim;
%for i=1:size(timpact,2)
%    impactline=[timpact(i) timpact(i)];
%    plot(impactline,[lim(1) lim(2)],'-.'); hold on
%end


%animate
