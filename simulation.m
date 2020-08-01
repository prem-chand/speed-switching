clear;
close all;
clc;
M=6;
te=0;
tplot=0;timpact=0;
ntouch=0;ptouch=0;
accel(1)=0;
simout.t=[];simout.stk=[];simout.h=[]; simout.swk=[]; simout.swf=[]; simout.head=[];
epsl=0.1;
tF=[0 1000];
Fext=[0   0]';
d = 0;
num_steps = 10;
P =[0.1313         0         0         0    0.0031         0         0         0
    0    0.1313         0         0         0    0.0031         0         0
    0         0    0.1313         0         0         0    0.0031         0
    0         0         0    0.1313         0         0         0    0.0031
    0.0031         0         0         0    0.0066         0         0         0
    0    0.0031         0         0         0    0.0066         0         0
    0         0    0.0031         0         0         0    0.0066         0
    0         0         0    0.0031         0         0         0    0.0066];
% P =[1.3125         0         0         0    0.0313         0         0         0;
%          0    1.3125         0         0         0    0.0313         0         0;
%          0         0    1.3125         0         0         0    0.0313         0;
%          0         0         0    1.3125         0         0         0    0.0313;
%     0.0313         0         0         0    0.0664         0         0         0;
%          0    0.0313         0         0         0    0.0664         0         0;
%          0         0    0.0313         0         0         0    0.0664         0;
%          0         0         0    0.0313         0         0         0    0.0664];

fix=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0_slow.mat');
% fix = open('optimization_red/fp1-2p5.mat');
% fix = open('optimization_red/fixed_point_guess1.mat');
% ival=fix.fpStates';
%fix=open('/Users/mohamadshafiee/Desktop/Reasearch/programming/5 link RABBIT/fixedpointforfivelink/fixed_stairs3.mat');
ival=fix.x;
% ival=fix.fpStates';

clear  theta_minus theta_plus betta XX event
global theta_minus theta_plus betta event M_bez
M_bez=6;
XX=[];
XX=ival(1:10)';

%event=XX(2);
%XX(:,1)=[2.91;3.54;-0.25;-0.41;0.056;-0.5;0.3;-0.3;-2.1;-0.56]    %book initial value before impact
%XX(:,1)=[3.54;2.91;-0.41;-0.25;0.056;0.3;-0.5;-2.1;-0.3;-0.56];     %book initial value after impact

c=[1 1 0 0.5 0];
[betta,theta_minus,theta_plus]=fcn_alpha_red(ival');
options = odeset('Events',@(t,q)touchdown5_slope(t,q,d),'RelTol',1e-5,'AbsTol',1e-4);
[t,y,te,ye,ie] = ode45(@(t,Z) return_map5_original(t,Z,Fext,tF,betta,theta_minus,theta_plus,M_bez),[0 20],XX,options);

zeta_fp = gamma0(y(end,:))*y(end,6:10)';
% XX = XX+rand(10,1)/50;
% % XX = [0.0245    2.7261    3.1019    0.3860    0.4058    2.0350   -1.7575   -2.0389    0.3192    0.8157]';
% XX = [0.0266    2.7227    3.0870    0.3662    0.4139    1.6380   -1.3666   -1.6606    0.1610    0.6865]';

% XX = [0.0313    3.1631    2.7778    0.2785    0.3782    2.8867   -2.3318   -1.5657    1.1593   -9.3207]';
% norm(XX-ival(1:10)')
% XX = [0.0368    3.1503    2.7712    0.2850    0.3269    2.7673   -2.2686   -0.7379    1.4193   -9.2453]';

%XX(1,1)=XX(1,1)-0.05;  %perturbation
%XX(:,1)=[0.1001;2.8435;3.063;0.175;0.17;3.103;-2.143;-3.072;0.154;1.96]
%XX(6:10,1)=[3.49354452125262;-2.41336795730199;-3.45828796801275;0.176033366955347;2.20745934941741;]
%XX(6:10,1)=[3.54502320297815;-2.44931414232983;-3.50908625540315;0.179478066319619;2.23975483281005;]; %new fixed point09 with the presence of force
%XX(6:10,1)=[3.56185549304056;-2.46069327873592;-3.52569850395199;0.179820021237419;2.25028180666301;];
%betta
%XX
%cost_verify=fixed_point_step([ival(1:4) ival(6:end)],Fext,tF)

% d = [-0.2863    0.3253   -0.4370   -0.5392    0.4223    0.2491    0.1812    0.3209   -0.9049   -0.3024]*1e-3;
d = d*ones(1,num_steps);
for i=1:num_steps
    % d = rand()*1e-3*1.0;
    % d = 0;
    clear s
    theta_minus=XX(1,i)+XX(3,i)+0.5*XX(5,i); %c*q-
    theta_plus=c*XX(1:5,i);                  %c*q+
    
    %thetadot_plus=c*XX(6:10,i)        %c*q.+
    %med=H0*XX(6:10,i)/thetadot_plus/M*(theta_minus-theta_plus)
    %betta(:,2)=med(1:4)+betta(:,1)   %equation 6.16a
    % Fext(:,i) = [-4 ; 0]+8*[rand() ; 0];
    Fext(:,i) = [0 ; 0];
    options = odeset('Events',@(t,q)touchdown5_slope(t,q,d(i)),'RelTol',1e-5,'AbsTol',1e-4);
    [t,y,te,ye,ie] = ode45(@(t,Z) return_map5_original(t+tplot,Z,Fext(:,i),tF,betta,theta_minus,theta_plus,M_bez),[0 20],XX(:,i),options);
    % Fext
    %XX(:,i)=ye'
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
    
    [XX(:,i+1) , Fimp] = impact_map5(ye);
    timpact(i)=tplot;
    % Fground(1)/Fground(2)
    
    %plotting states
    impactline=[timpact(i) timpact(i)];
    % figure(1)
    % for j=1:10
    %     subplot(2,5,j); plot(t+tplot,y(:,j),'LineWidth',1);xlabel('Time (sec)','fontsize',12);ylabel(['q',num2str(j)],'fontsize',14);%title(['q',num2str(j)]);
    %     hold on; plot(impactline,[min(y(:,j)) max(y(:,j))],'-.'); hold on
    % end
    
    %plotting q versus qdot
    %figure(8)
    % plot(y(:,8),y(:,3));xlabel('q1^.','fontsize',12);ylabel('q1','fontsize',12);hold on
    
    %plotting theta
    % figure(2)
    %     subplot(2,1,1); plot(t'+tplot,c*y(:,1:5)');xlabel('Time (sec)','fontsize',12);ylabel('\theta','fontsize',14);lim=ylim;
    %     %hold on; plot(impactline,[min(c*y(:,1:5)') max(c*y(:,1:5)')],'-.'); hold on
    %     hold on; plot(impactline,[lim(1) lim(2)],'-.'); hold on
    %     subplot(2,1,2); plot(t'+tplot,c*y(:,6:10)');xlabel('Time (sec)','fontsize',12);ylabel('\theta^.','fontsize',14);lim=ylim;
    %     hold on; plot(impactline,[lim(1) lim(2)],'-.'); hold on
    
    %plotting theta versus thetadot
    % figure(3)
    %     plot(c*y(:,1:5)',c*y(:,6:10)');xlabel('\theta','fontsize',14);ylabel('\theta^.','fontsize',14);hold on
    %     %plot([c*y(1,1:5)',c*y(size(y,1),1:5)'],[c*y(1,6:10)',c*y(size(y,1),6:10)']);
    %     plot([c*XX(1:5,i+1),c*y(size(y,1),1:5)'],[c*XX(6:10,i+1),c*y(size(y,1),6:10)']);
    
    %computing inputs outputs and nomralized theta (s)
    clear u h hdot GRF s;
    for tp=1:size(t)
        [fx,gx,gforce]=fcn_stance_dynamics5(y(tp,:));
        % Fext = force(t(tp)+tplot);
        % Fext_plot = [interp1(tF,Fext,t(tp)+tplot);0];
        u(tp,:)=fcn_stance_controller5(y(tp,:),fx,gx,gforce,[0;0],betta,theta_minus,theta_plus,M_bez);
        if M_bez==6
            [h(tp,:),hdot(tp,:)]=fcn_liederivative(y(tp,:),betta,theta_minus,theta_plus,M_bez);
        elseif M_bez==5
            [h(tp,:),hdot(tp,:)]=fcn_liederivative_m5(y(tp,:),betta);
        end
        s(tp)=(c*y(tp,1:5)'-theta_plus)/(theta_minus-theta_plus);
        % GRF(tp,:)=[cosd(angle) -sind(angle) ; sind(angle) cosd(angle)]*fcn_out_force(y(tp,:)',u(tp,:)');
        % friction_cone(tp) = abs(GRF(tp,1)/GRF(tp,2));
        cm(tp,:)=fcn_cm_kinematics(y(tp,1:5));
    end
    section_state(:,i) = y(end,:)';
    V_eta(i) = [(1/epsl)*h(end,:),hdot(end,:)]*P*[(1/epsl)*h(end,:),hdot(end,:)]';
    zeta = gamma0(y(end,:))*y(end,6:10)';
    V_xi(i) = (zeta-zeta_fp)^2;
    V(i) = 100*V_eta(i) + V_xi(i)
    output_error=[h hdot];
    norm_output_error = max(sqrt(sum(output_error.^2,2)));
    % GRF(tp,:)=fcn_GRF(y(tp,:),u(tp,:));
    % Fimp = [cosd(angle) -sind(angle) ; sind(angle) cosd(angle)]*Fimp';
    % friction_cone = [friction_cone abs(Fimp(1)/Fimp(2))];
    % coeff_fric = max(friction_cone);
    
    % plotting inputs
    % figure(4)
    % subplot(2,2,1);plot(t+tplot,u(:,1));ylabel('u_1 (Nm)','fontsize',16);
    % hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    % subplot(2,2,2);plot(t+tplot,u(:,2));ylabel('u_2 (Nm)','fontsize',16);
    % hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    % subplot(2,2,3);plot(t+tplot,u(:,3));xlabel('Time (s)','fontsize',16);ylabel('u_3 (Nm)','fontsize',16);
    % hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    % subplot(2,2,4);plot(t+tplot,u(:,4));xlabel('Time (s)','fontsize',16);ylabel('u_4 (Nm)','fontsize',16);
    % hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    u_max(i) = max(max(abs(u)))
    
    %plotting outputs
    figure(5)
    subplot(2,4,1);plot(t+tplot,h(:,1));title('h1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    subplot(2,4,2);plot(t+tplot,h(:,2));title('h2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    subplot(2,4,3);plot(t+tplot,h(:,3));title('h3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    subplot(2,4,4);plot(t+tplot,h(:,4));title('h4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    
    subplot(2,4,5);plot(t+tplot,hdot(:,1));title('h^.1');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    subplot(2,4,6);plot(t+tplot,hdot(:,2));title('h^.2');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    subplot(2,4,7);plot(t+tplot,hdot(:,3));title('h^.3');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    subplot(2,4,8);plot(t+tplot,hdot(:,4));title('h^.4');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
    %plotting s
    %figure(4)
    %plot(t+tplot,s);title('Normalized Theta');hold on;plot(impactline,[min(s) max(s)],'-.'); hold on
    
    %plotting velocity of hip
    %figure(6)
    %plot(t+tplot,vh);title('Velocity of Hip');hold on;lim=ylim;plot(impactline,[lim(1) lim(2)],'-.'); hold on
    
    % figure(7)
    % plot(t+tplot,Menergy);title('Mechanical Energy');hold on;lim=ylim;plot(impactline,[lim(1) lim(2)],'-.'); hold on
    %
    % energy_loss(i)=Menergy(size(Menergy,1))-fcn_energy(XX(:,i+1));
    % figure(8)
    % plot(i,energy_loss(i),'*');title('Energy Loss');hold on
    
    step_length(i,:)=fcn_position_swingfoot(ye);
    average_speed(i)=step_length(i,1)/te
    
    if i>=2
        accel(i)=(average_speed(i)-average_speed(i-1))/te;
    end
    tplot=tplot+te;
    
    error=norm(XX(:,i+1)-XX(:,i));
    % end
    
    %thetadot_minus=c*ye(6:10)';
    %checkpoint= betta(:,7)- H0*ye(6:10)'*(theta_minus-theta_plus)/thetadot_minus/M
    %betta
    
end
XX;
timpact(i+1)=tplot;
hold on;
figure(1)
hold on;
plot([1:num_steps],V,'x','color','r','linewidth',1.5);
% xlabel('Step Number','fontsize',30);
% ylabel('V_x','fontsize',30);
figure(2)
plot([1:num_steps],u_max,'x','linewidth',2);
hold on;
% plot([1:num_steps],20*ones(1,num_steps),'--','color','k','linewidth',1);
hold off;

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
