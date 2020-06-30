clear;
close all;
clc;
M=6;
c_plus=[1 1 0 0.5 0];
c_minus=[1 0 1 0 0.5];
ifenter=1;
te=0;
tplot=0;timpact=0;
ntouch=0;ptouch=0;
simout.t=[];simout.stk=[];simout.h=[]; simout.swk=[]; simout.swf=[]; simout.head=[];
accel(1)=0;
F_step=0;
c=[1 1 0 0.5 0];

%t0f=0.3;t1f=3;
%tF=[0  1 1.0001 8 8.001 20];
%Fext=[0 0 50    50 0       0];
%tF=  [0  t0f t0f+0.0001 t1f t1f+0.0001 30];
%Fext=[0   0    10        10      0   0];

F_index=1;
ctr_index=3;
tI=0.2241;% 0.6357
%tI=0.5651;
num_steps=10;

%tF=  [0 tI*F_index tI*F_index+0.01 50];
%Fext=[0     0          0         0];
%tF=zeros(1,2*(num_steps+1))
%Fext=zeros(1,2*(num_steps+1))
tF=[0 .2241 .2242 0.4379 0.4380  50];
Fext=[0 0    10   10     0      0];
%tF=[0 50];
%Fext=[0  0];
global XX_noforce

fix=open('/Users/mohamadshafiee/Desktop/Reasearch/programming/5 link RABBIT/fixedpointforfivelink/fixed_point_9.mat');
ival_1=fix.fpStates
%fix=open('/Users/mohamadshafiee/Desktop/Reasearch/programming/5 link RABBIT/fixedpointforfivelink/fixed2_10N.mat');
%ival_1=fix.x;
XX1=ival_1(1:10)
XX2=XX_noforce(:,1)

fix=open('/Users/mohamadshafiee/Desktop/Reasearch/programming/5 link RABBIT/fixedpointforfivelink/converg_5.mat');
ival_2=[XX2' fix.x]'



[betta1,theta_minus_1,theta_plus_1]=fcn_alpha_red(ival_1)
[betta2,theta_minus_2,theta_plus_2]=fcn_alpha_red(ival_2)

clear  theta_minus theta_plus betta event XX
global theta_minus theta_plus betta event

XX=XX1;
betta=betta1;
theta_minus=theta_minus_1;
theta_plus=theta_plus_1;   
event=XX1(2);

for i=1:num_steps
clear s

if i==ctr_index+1
    betta=betta2;
    %[betta,theta_minus,theta_plus,event]=fcn_transition_ctr(ival_1,ival_2);

elseif i==ctr_index+2 
    betta=betta1;
    
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
ptouch=ptouch+pswf(size(pswf,1),1);
str_lng(i,1)=ntouch;
str_lng(i,2)=ptouch;

XX(:,i+1)=impact_map5(ye);
timpact(i)=tplot;

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
    plot(c*y(:,1:5)',c*y(:,6:10)');xlabel('\theta','fontsize',14);ylabel('\theta^.','fontsize',14);hold on
    %plot([c*y(1,1:5)',c*y(size(y,1),1:5)'],[c*y(1,6:10)',c*y(size(y,1),6:10)']);
    plot([c*XX(1:5,i+1),c*y(size(y,1),1:5)'],[c*XX(6:10,i+1),c*y(size(y,1),6:10)']);
%computing inputs outputs and nomralized theta (s)
clear u h hdot Ftoe;
for tp=1:size(t)
[fx,gx,gforce]=fcn_stance_dynamics5(y(tp,:));
Fext_plot = [interp1(tF,Fext,t(tp));0];
u(tp,:)=fcn_stance_controller5(y(tp,:),fx,gx,gforce,Fext_plot);
Ftoe(tp,:)=fcn_out_force(y(tp,:),u(tp,:));
[h(tp,:),hdot(tp,:)]=fcn_liederivative(y(tp,:),betta);
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
subplot(2,4,1);plot(t+tplot,h(:,1));ylabel('y_1(rad)','fontsize',16);hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
subplot(2,4,2);plot(t+tplot,h(:,2));ylabel('y_2(rad)','fontsize',16);hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
subplot(2,4,3);plot(t+tplot,h(:,3));ylabel('y_3(rad)','fontsize',16);hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.'); 
subplot(2,4,4);plot(t+tplot,h(:,4));ylabel('y_4(rad)','fontsize',16);hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');

subplot(2,4,5);plot(t+tplot,hdot(:,1));xlabel('Time (s)','fontsize',16);ylabel('y^._1(rad/s)','fontsize',16);hold on;lim=ylim; plot(impactline,[lim(1) lim(2)],'-.');
subplot(2,4,6);plot(t+tplot,hdot(:,2));xlabel('Time (s)','fontsize',16);ylabel('y^._2(rad/s)','fontsize',16);hold on; plot(impactline,[lim(1) lim(2)],'-.');
subplot(2,4,7);plot(t+tplot,hdot(:,3));xlabel('Time (s)','fontsize',16);ylabel('y^._3(rad/s)','fontsize',16);hold on; plot(impactline,[lim(1) lim(2)],'-.'); 
subplot(2,4,8);plot(t+tplot,hdot(:,4));xlabel('Time (s)','fontsize',16);ylabel('y^._4(rad/s)','fontsize',16);hold on; plot(impactline,[lim(1) lim(2)],'-.');
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
average_speed(i)=step_length(i,1)/te;


if i>=2
    accel(i)=(average_speed(i)-average_speed(i-1))/te;
end
tplot=tplot+te;


end
XX
timpact(i+1)=tplot;

figure (8)
subplot(1,1,1);plot(timpact(2:end),average_speed,'rs');xlim([0 tplot]);xlabel('Time (s)','fontsize',16);ylabel('Average Speed (m/s)','fontsize',16);
%subplot(3,1,2);plot(timpact(2:end),accel,'rs');xlim([0 tplot]);title('Average Accelaration');hold on; 
%subplot(3,1,3);plot(tF,Fext,'--rs');xlim([0 tplot]);title('External Force');
F_mean=mean(Fext)

mean_accel=(average_speed(i)-average_speed(1))/(tplot-timpact(2))

norm(XX(:,1)-XX(:,5))
%animate
