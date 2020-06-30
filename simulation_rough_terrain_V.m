clear;
% close all;
clc;
M=6;
te=0;
tplot=0;timpact=0;
ntouch=0;ptouch=0;
accel(1)=0;
simout.t=[];simout.stk=[];simout.h=[]; simout.swk=[]; simout.swf=[]; simout.head=[];
tF=[0 500];
Fext=[0   0];
angle = 0;
num_steps = 20;
P =[1.3125         0         0         0    0.0313         0         0         0;
         0    1.3125         0         0         0    0.0313         0         0;
         0         0    1.3125         0         0         0    0.0313         0;
         0         0         0    1.3125         0         0         0    0.0313;
    0.0313         0         0         0    0.0664         0         0         0;
         0    0.0313         0         0         0    0.0664         0         0;
         0         0    0.0313         0         0         0    0.0664         0;
         0         0         0    0.0313         0         0         0    0.0664];
    
fix=open('fixedpointforfivelink/fixed_point_0deg_1_cone.mat');
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
options = odeset('Events',@(t,q)touchdown5(t,q,angle),'RelTol',1e-5,'AbsTol',1e-4);
[t,y,te,ye,ie] = ode45(@(t,Z) return_map5(t,Z,Fext,tF,betta,theta_minus,theta_plus,M_bez),[0 20],XX,options);

zeta_fp = gamma0(y(end,:))*y(end,6:10)';
% XX = XX+rand(10,1)/100;
XX = [0.0313    3.1631    2.7778    0.2785    0.3782    2.8867   -2.3318   -1.5657    1.1593   -9.3207]';
norm(XX-ival(1:10)')
% XX = [0.0368    3.1503    2.7712    0.2850    0.3269    2.7673   -2.2686   -0.7379    1.4193   -9.2453]';

%XX(1,1)=XX(1,1)-0.05;  %perturbation
%XX(:,1)=[0.1001;2.8435;3.063;0.175;0.17;3.103;-2.143;-3.072;0.154;1.96]
%XX(6:10,1)=[3.49354452125262;-2.41336795730199;-3.45828796801275;0.176033366955347;2.20745934941741;]
%XX(6:10,1)=[3.54502320297815;-2.44931414232983;-3.50908625540315;0.179478066319619;2.23975483281005;]; %new fixed point09 with the presence of force
%XX(6:10,1)=[3.56185549304056;-2.46069327873592;-3.52569850395199;0.179820021237419;2.25028180666301;];
%betta
%XX
%cost_verify=fixed_point_step([ival(1:4) ival(6:end)],Fext,tF)

for i=1:num_steps
clear s
theta_minus=XX(1,i)+XX(3,i)+0.5*XX(5,i); %c*q-
theta_plus=c*XX(1:5,i);                  %c*q+

% Fext(:,i) = [1*rand() ; 0];
% Fext = [0 ; 0];
options = odeset('Events',@(t,q)touchdown5(t,q,angle),'RelTol',1e-5,'AbsTol',1e-4);
[t,y,te,ye,ie] = ode45(@(t,Z) return_map5(t+tplot,Z,Fext,tF,betta,theta_minus,theta_plus,M_bez),[0 20],XX(:,i),options);
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
impactline=[timpact(i) timpact(i)];

clear u h hdot GRF s;
for tp=1:size(t)
[fx,gx,gforce]=fcn_stance_dynamics5(y(tp,:));
Fext = force(t(tp)+tplot);
% Fext_plot = [interp1(tF,Fext,t(tp)+tplot);0];
% u(tp,:)=fcn_stance_controller5(y(tp,:),fx,gx,gforce,Fext_plot,betta,theta_minus,theta_plus,M_bez);
if M_bez==6
    [h(tp,:),hdot(tp,:)]=fcn_liederivative(y(tp,:),betta,theta_minus,theta_plus,M_bez);
elseif M_bez==5
    [h(tp,:),hdot(tp,:)]=fcn_liederivative_m5(y(tp,:),betta);
end
s(tp)=(c*y(tp,1:5)'-theta_plus)/(theta_minus-theta_plus);

cm(tp,:)=fcn_cm_kinematics(y(tp,1:5));
end
section_state(:,i) = y(end,:)';
V_eta(i) = [h(end,:),hdot(end,:)]*P*[h(end,:),hdot(end,:)]'
zeta = gamma0(y(end,:))*y(end,6:10)';
V_xi(i) = (zeta-zeta_fp)^2
V(i) = 10*V_eta(i) + V_xi(i);
output_error=[h hdot];
norm_output_error = max(sqrt(sum(output_error.^2,2)));
% GRF(tp,:)=fcn_GRF(y(tp,:),u(tp,:));
Fimp = [cosd(angle) -sind(angle) ; sind(angle) cosd(angle)]*Fimp';

step_length(i,:)=fcn_position_swingfoot(ye);
average_speed(i)=step_length(i,1)/te;

if i>=2
    accel(i)=(average_speed(i)-average_speed(i-1))/te;
end
tplot=tplot+te;

    error=norm(XX(:,i+1)-XX(:,i))


end
XX;
timpact(i+1)=tplot;
figure(2)
plot([1:num_steps],V_xi,'x','linewidth',2);
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
