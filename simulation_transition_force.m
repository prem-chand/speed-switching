clear;
close all;
clc;
M=6;
c_plus=[1 1 0 0.5 0];
c_minus=[1 0 1 0 0.5];
ifenter=1;
te=0;
tplot=0;timpact=0;
ntouch=0;ptouch=[0 0];
simout.t=[];simout.stk=[];simout.h=[]; simout.swk=[]; simout.swf=[]; simout.head=[];
accel(1)=0;
F_step=0;

%t0f=0.3;t1f=3;
%tF=[0  1 1.0001 8 8.001 20];
%Fext=[0 0 50    50 0       0];
%tF=  [0  t0f t0f+0.0001 t1f t1f+0.0001 30];
%Fext=[0   0    10        10      0   0];

ctr_index=[10 90];
tI=0.2241;% 0.6357
%tI=0.5651;
num_steps=3;
a_des=[0 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];  % must have dimension=num_steps
%a_des=[0 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0 0 0 0 0 0 0];
%a_des=zeros(1,num_steps);
%Pu=0.5;
%kpid=[110*0.6 1.2*55/Pu 0.6/8*55*Pu*0];

kpid=[80 0 0];
tF=[0 50];
Fext=[0  0];

%tF=  [0 tI*F_index tI*F_index+0.01 50];
%Fext=[0     0          0         0];
%tF=zeros(1,2*(num_steps+1))
%Fext=zeros(1,2*(num_steps+1))

%tF=  [0 0.2241 0.2242  10 ];
%Fext=[0   0       10   10 ];



%fix=open('/Users/mohamadshafiee/Desktop/Reasearch/programming/5 link RABBIT/fixedpointforfivelink/fixed_point_9.mat');
%ival_1=fix.fpStates;
fix=open('fixedpointforfivelink/fixed_point_96.mat');
ival_1=fix.x';
fix=open('fixedpointforfivelink/fixed_point_2.mat');
ival_2=fix.fpStates;

XX1=ival_1(1:10);
XX2=ival_2(1:10);
clear  theta_minus theta_plus betta event XX M_bez
global theta_minus theta_plus betta event M_bez

[betta1,theta_minus_1,theta_plus_1]=fcn_alpha_red(ival_1)
[betta2,theta_minus_2,theta_plus_2]=fcn_alpha_red(ival_2)
M1=(max(size(ival_1))-10)/4+2;
M2=(max(size(ival_2))-10)/4+2;
M_bez=M1;

XX=XX1;
betta=betta1;
theta_minus=theta_minus_1;
theta_plus=theta_plus_1;   
event=XX1(2);

for i=1:num_steps
clear s

if i==ctr_index(1)+1
    [betta,theta_minus,theta_plus,event]=fcn_transition_ctr(ival_1,ival_2);
    M_bez=6;
    
elseif (i>ctr_index(1)+1) & (ifenter==1)
    theta_minus = theta_minus_2;
    theta_plus  = theta_plus_2;
    betta=betta2;
    M_bez=M2;
    ifenter=0;
end
if i==ctr_index(2)+1
    [betta,theta_minus,theta_plus,event]=fcn_transition_ctr(ival_2,ival_1);
    M_bez=6;
    
elseif (i>ctr_index(2)+1) & (ifenter==0)
    theta_minus = theta_minus_1;
    theta_plus  = theta_plus_1;
    betta=betta1;
    M_bez=M1;
    ifenter=1;
end
     
    


%thetadot_plus=c*XX(6:10,i)        %c*q.+
%med=H0*XX(6:10,i)/thetadot_plus/M*(theta_minus-theta_plus)
%betta(:,2)=med(1:4)+betta(:,1)   %equation 6.16a

options = odeset('Events',@touchdown5,'RelTol',1e-5,'AbsTol',1e-5);
[t,y,te,ye,ie] = ode45(@(t,Z) return_map5(t+tplot,Z,Fext,tF),[0 20],XX(:,i),options);
%XX(:,i)=ye'

simout.t=[simout.t;t+tplot];
[pstk ph pswk pswf phead vh]=out_kinematics(y);
%simout.stk = [simout.stk  ; pstk(:,1)+str_lng pstk(:,2)];
%simout.h   = [simout.h    ; ph(:,1)+str_lng ph(:,2)];
%simout.swk = [simout.swk  ; pswk(:,1)+str_lng pswk(:,2)];
%simout.swf = [simout.swf  ; pswf(:,1)+str_lng pswf(:,2)];
%simout.head= [simout.head ; phead(:,1)+str_lng phead(:,2)];
simout.stk = [simout.stk  ; pstk];
simout.h   = [simout.h    ; ph];
simout.swk = [simout.swk  ; pswk];
simout.swf = [simout.swf  ; pswf];
simout.head= [simout.head ; phead];

ntouch=ntouch+size(pswf,1);
%ptouch=ptouch+pswf(size(pswf,1));
ptouch=ptouch+pswf(end,:);

str_lng(i,1)=ntouch;
str_lng(i,2:3)=ptouch;

XX(:,i+1)=impact_map5(ye);
timpact(i)=tplot;
tstep(i)=te;

%plotting states
impactline=[timpact(i) timpact(i)];
figure(1)
for j=1:10
    subplot(2,5,j); plot(t+tplot,y(:,j),'LineWidth',1);xlabel('Time (sec)','fontsize',12);ylabel(['q',num2str(j)],'fontsize',14);%title(['q',num2str(j)]);
    hold on; plot(impactline,[min(y(:,j)) max(y(:,j))],'-.'); hold on
end

%plotting q versus qdot
%figure(8)
% plot(y(:,8),y(:,3));xlabel('q1^.','fontsize',12);ylabel('q1','fontsize',12);hold on

%plotting theta
figure(2)
    subplot(2,1,1); plot(t'+tplot,c_plus*y(:,1:5)');xlabel('Time (sec)','fontsize',12);ylabel('\theta','fontsize',14);lim=ylim;
    %hold on; plot(impactline,[min(c*y(:,1:5)') max(c*y(:,1:5)')],'-.'); hold on
    hold on; plot(impactline,[lim(1) lim(2)],'-.'); hold on
    subplot(2,1,2); plot(t'+tplot,c_plus*y(:,6:10)');xlabel('Time (sec)','fontsize',12);ylabel('\theta^.','fontsize',14);lim=ylim;
    hold on; plot(impactline,[lim(1) lim(2)],'-.'); hold on
    
%plotting theta versus thetadot
figure(3)
    plot(c_plus*y(:,1:5)',c_plus*y(:,6:10)');xlabel('\theta','fontsize',14);ylabel('\theta^.','fontsize',14);hold on
    %plot([c_plus*y(1,1:5)',c_plus*y(size(y,1),1:5)'],[c_plus*y(1,6:10)',c_plus*y(size(y,1),6:10)']);
     plot([c_plus*y(size(y,1),1:5)',c_plus*XX(1:5,i+1)],[c_plus*y(size(y,1),6:10)',c_plus*XX(6:10,i+1)]);
%computing inputs outputs and nomralized theta (s)
clear u h hdot Ftoe;
for tp=1:size(t)
[fx,gx,gforce]=fcn_stance_dynamics5(y(tp,:));
%t(tp)+tplot
Fext_plot = [interp1(tF,Fext,t(tp)+tplot);0];
u(tp,:)=fcn_stance_controller5(y(tp,:),fx,gx,gforce,Fext_plot);
Ftoe(tp,:)=fcn_out_force(y(tp,:),u(tp,:));
if M_bez==6
    [h(tp,:),hdot(tp,:)]=fcn_liederivative(y(tp,:),betta);
elseif M_bez==5
    [h(tp,:),hdot(tp,:)]=fcn_liederivative_m5(y(tp,:),betta);
end
s(tp)=(c_plus*y(tp,1:5)'-theta_plus)/(theta_minus-theta_plus);
end

%plotting inputs
figure(4)
subplot(2,2,1);plot(t+tplot,u(:,1));ylabel('u_1 (Nm)','fontsize',16);
hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
subplot(2,2,2);plot(t+tplot,u(:,2));ylabel('u_2 (Nm)','fontsize',16);
hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
subplot(2,2,3);plot(t+tplot,u(:,3));xlabel('Time (s)','fontsize',16);ylabel('u_3 (Nm)','fontsize',16);
hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.'); 
subplot(2,2,4);plot(t+tplot,u(:,4));xlabel('Time (s)','fontsize',16);ylabel('u_4 (Nm)','fontsize',16);
hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.'); 

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
figure(6)
plot(t+tplot,s);title('Normalized Theta');hold on;plot(impactline,[min(s) max(s)],'-.'); hold on

%plotting velocity of hip
%figure(7)
%plot(t+tplot,vh);title('Velocity of Hip');hold on;lim=ylim;plot(impactline,[lim(1) lim(2)],'-.'); hold on

%plot reaction forces
figure(7)
subplot(3,1,1);plot(t+tplot,Ftoe(:,1));title('FT');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
subplot(3,1,2);plot(t+tplot,Ftoe(:,2));title('FN');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
subplot(3,1,3);plot(t+tplot,Ftoe(:,1)/Ftoe(:,2));title('\mu');hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');


step_length(i,:)=fcn_position_swingfoot(ye);
robot_speed(i)=step_length(i,1)/te;
robot_pos(i)=sum(step_length(:,1));
leader_speed(1)=robot_speed(1);
leader_pos(1)=robot_pos(1);
if i>=2
    accel(i)=(robot_speed(i)-robot_speed(i-1))/te;
    leader_speed(i)=leader_speed(i-1)+a_des(i)*te;
    leader_pos(i)=leader_pos(i-1)+0.5*a_des(i)*te^2+leader_speed(i-1)*te;

end
tplot=tplot+te;

%if (i>=ctr_index+3) | (i<ctr_index)
%    F_step=35*(a_des(i)-accel(i))+F_step;
%end

%force control through accelaration
error(i)=a_des(i)-accel(i);
if i==1 
    F_step=kpid(1)*error(i);
%elseif (i>=ctr_index+3) | (i<ctr_index)
else
    F_step=kpid(1)*error(i)+kpid(2)*sum(error)+kpid(3)*(accel(i)-accel(i-1))+F_step
end
%end
%force control through position
%if i==1 
%    F_step=kpd(1)*(leader_pos(i)-robot_pos(i))+F_step;
%elseif (i>=ctr_index+3) | (i<ctr_index)
%    F_step=kpd(1)*(leader_pos(i)-robot_pos(i))+kpd(2)*(leader_pos(i)-leader_pos(i-1)-robot_pos(i)+robot_pos(i-1))+F_step;
%end

%if a_des(i)==0 
%    F_step=0;
%end

% step-like force
%Fext(2*i+1)=F_step;Fext(2*i+2)=F_step;
%tF(2*i)=tplot;tF(2*i+1)=tplot+0.0001;tF(2*i+2)=tplot+3*te;
%
% slope_like force
%Fext(i+2)=F_step;
%tF(i+1)=tplot;tF(i+2)=tplot+1.2*te;
H0=[zeros(4,1) eye(4)];
thetadot_minus=c_plus*ye(6:10)';
checkpoint(:,i)= betta(:,7)- H0*ye(6:10)'*(theta_minus-theta_plus)/thetadot_minus/M_bez;   


end
XX
timpact(i+1)=tplot
step_length
tstep
%figure (8)
%subplot(4,1,1);plot(tF,Fext,'-rs','LineWidth',2);xlim([0 tplot]);ylabel('External Force (N)','fontsize',12);
%subplot(4,1,2);plot(timpact(2:end),robot_pos*100,'b*',timpact(2:end),leader_pos*100,'rs');xlim([0 tplot]);ylabel('Position (cm)','fontsize',12);legend('Robot','Leader');
subplot(1,1,1);plot(timpact(2:end),robot_speed,'b*',timpact(2:end),leader_speed,'rs');xlim([0 tplot]);ylabel('Speed (m/s)','fontsize',18);legend('Robot','Leader');xlabel('Time (s)','fontsize',18);
%subplot(4,1,4);plot(timpact(2:end),accel,'b*');xlim([0 tplot]);hold on; plot(timpact(2:end),a_des,'rs');hold on;xlabel('Time (s)','fontsize',16);ylabel('Accelaration (m/s^2)','fontsize',12);legend('Robot','Leader');
%F_mean=mean(Fext)

mean_accel=(robot_speed(i)-robot_speed(1))/(tplot-timpact(2))
%sum(step_length(2:end,1));
%mean_accel=(sum(step_length(2:end,1))-robot_speed(1)*(tplot-timpact(2)))*2/(tplot-timpact(2))^2
%cost=norm(a_des-accel)
checkpoint
betta(:,6)

%animate
