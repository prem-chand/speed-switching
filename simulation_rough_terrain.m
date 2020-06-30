clear;
close all;
clc;
M=6;
te=0;
tplot=0;timpact=0;
ntouch=0;ptouch=0;
accel(1)=0;
simout.t=[];simout.stk=[];simout.h=[]; simout.swk=[]; simout.swf=[]; simout.head=[];
%[4.60950011028077,2.69622282334424,2.99888042275339,4.12534914292605,1.15725818250891;4.35950011028077,3.07122282334424,2.99888042275339,4.12534914292605,0.907258182508908;4.35950011028077,3.13372282334424,2.99888042275339,4.12534914292605,0.907258182508908;4.35950011028077,3.07122282334424,2.99888042275339,4.12534914292605,0.907258182508908;4.60950011028077,3.19622282334424,2.99888042275339,4.12534914292605,1.15725818250891;4.10950011028077,3.13372282334424,2.99888042275339,4.12534914292605,0.907258182508908;4.35950011028077,3.13372282334424,2.99888042275339,4.12534914292605,0.407258182508908;4.60950011028077,2.69622282334424,2.99888042275339,4.12534914292605,1.15725818250891;4.35950011028077,3.13372282334424,2.99888042275339,4.12534914292605,1.15725818250891;4.60950011028077,1.69622282334424,2.99888042275339,5.12534914292605,2.15725818250891;]
%t0f=0.3;t1f=3;
%tF=[0  0.2241*3 0.2241*3+0.001 0.2241*5 0.2241*5+0.001 20];
%Fext=[0  0        10             10       0      0];
%tF=  [0  t0f t0f+0.0001 t1f t1f+0.0001 30];
%Fext=[0   0    10        10      0   0];
%
%tF=[0 .2241 .2242 0.4379 0.4380  50];
%Fext=[0 0    10   10     0      0];
tF=[0 50];
Fext=[0   0];
%a_des=[0 0.1 0.1 0.1 0.1 0.1 0.1];
%optimized only for force
%tF=[0,0.224051700858414,0.448114997580434,0.672186246208195,0.896265325583177,1.12035190722882,1.34444580807905];
%Fext=[0 0 4.5 2.4 3.5 3.3 3.4];
%tF=[0,0.224051700858414,0.448114997580434,0.672186246208195,0.896265325583177,1.12035190722882,1.34444580807905,1.65];
%Fext=[0 0 0 4.6095    2.6962    2.9989    4.1253    1.1573];
%optimized for both Fext and tF
%par=[4.06254276640266,1.93976238300905,1.49920985312643,1.27202617444956,0.554558673432953,0.431209240089772,1.08157713051977,0.618192641635759,0.769877218774926,1.02870594212796];
%Fext=[0 0 0 par(1:5)]
%tF=[0 0.224 0.448 0.448+par(6) 0.448+par(6)+par(7) 0.448+par(6)+par(7)+par(8) 0.448+par(6)+par(7)+par(8)+par(9) 0.448+par(6)+par(7)+par(8)+par(9)+par(10)]

    
fix=open('fixedpointforfivelink/fixed_point_9_s.mat');
% ival=fix.fpStates';
%fix=open('/Users/mohamadshafiee/Desktop/Reasearch/programming/5 link RABBIT/fixedpointforfivelink/fixed_stairs3.mat');
ival=fix.x;
% ival=fix.fpStates';

clear  theta_minus theta_plus betta XX event
global theta_minus theta_plus betta event M_bez
M_bez=6;
XX=[];
XX=ival(1:10)';
event=XX(2);
%XX(:,1)=[2.91;3.54;-0.25;-0.41;0.056;-0.5;0.3;-0.3;-2.1;-0.56]    %book initial value before impact
%XX(:,1)=[3.54;2.91;-0.41;-0.25;0.056;0.3;-0.5;-2.1;-0.3;-0.56];     %book initial value after impact

c=[1 1 0 0.5 0];
[betta,theta_minus,theta_plus]=fcn_alpha_red(ival')


%XX(1,1)=XX(1,1)-0.05;  %perturbation
%XX(:,1)=[0.1001;2.8435;3.063;0.175;0.17;3.103;-2.143;-3.072;0.154;1.96]
%XX(6:10,1)=[3.49354452125262;-2.41336795730199;-3.45828796801275;0.176033366955347;2.20745934941741;]
%XX(6:10,1)=[3.54502320297815;-2.44931414232983;-3.50908625540315;0.179478066319619;2.23975483281005;]; %new fixed point09 with the presence of force
%XX(6:10,1)=[3.56185549304056;-2.46069327873592;-3.52569850395199;0.179820021237419;2.25028180666301;];
%betta
%XX
%cost_verify=fixed_point_step([ival(1:4) ival(6:end)],Fext,tF)

for i=1:5
clear s
theta_minus=XX(1,i)+XX(3,i)+0.5*XX(5,i); %c*q-
theta_plus=c*XX(1:5,i);                  %c*q+

%thetadot_plus=c*XX(6:10,i)        %c*q.+
%med=H0*XX(6:10,i)/thetadot_plus/M*(theta_minus-theta_plus)
%betta(:,2)=med(1:4)+betta(:,1)   %equation 6.16a

options = odeset('Events',@touchdown5,'RelTol',1e-4,'AbsTol',1e-3);
[t,y,te,ye,ie] = ode45(@(t,Z) return_map5(t+tplot,Z,Fext,tF),[0 20],XX(:,i),options);
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

XX(:,i+1)=impact_map5(ye);
timpact(i)=tplot

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
    subplot(2,1,1); plot(t'+tplot,c*y(:,1:5)');xlabel('Time (sec)','fontsize',12);ylabel('\theta','fontsize',14);lim=ylim;
    %hold on; plot(impactline,[min(c*y(:,1:5)') max(c*y(:,1:5)')],'-.'); hold on
    hold on; plot(impactline,[lim(1) lim(2)],'-.'); hold on
    subplot(2,1,2); plot(t'+tplot,c*y(:,6:10)');xlabel('Time (sec)','fontsize',12);ylabel('\theta^.','fontsize',14);lim=ylim;
    hold on; plot(impactline,[lim(1) lim(2)],'-.'); hold on
    
%plotting theta versus thetadot
figure(3)
    plot(c*y(:,1:5)',c*y(:,6:10)');xlabel('\theta','fontsize',14);ylabel('\theta^.','fontsize',14);hold on
    %plot([c*y(1,1:5)',c*y(size(y,1),1:5)'],[c*y(1,6:10)',c*y(size(y,1),6:10)']);
    plot([c*XX(1:5,i+1),c*y(size(y,1),1:5)'],[c*XX(6:10,i+1),c*y(size(y,1),6:10)']);

%computing inputs outputs and nomralized theta (s)
clear u h hdot GRF s;
for tp=1:size(t)
[fx,gx,gforce]=fcn_stance_dynamics5(y(tp,:));
Fext_plot = [interp1(tF,Fext,t(tp)+tplot);0];
u(tp,:)=fcn_stance_controller5(y(tp,:),fx,gx,gforce,Fext_plot);
if M_bez==6
    [h(tp,:),hdot(tp,:)]=fcn_liederivative(y(tp,:),betta);
elseif M_bez==5
    [h(tp,:),hdot(tp,:)]=fcn_liederivative_m5(y(tp,:),betta);
end
s(tp)=(c*y(tp,1:5)'-theta_plus)/(theta_minus-theta_plus);
GRF(tp,:)=fcn_out_force(y(tp,:)',u(tp,:)');
cm(tp,:)=fcn_cm_kinematics(y(tp,1:5));
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
%figure(4)
%plot(t+tplot,s);title('Normalized Theta');hold on;plot(impactline,[min(s) max(s)],'-.'); hold on

%plotting velocity of hip
%figure(6)
%plot(t+tplot,vh);title('Velocity of Hip');hold on;lim=ylim;plot(impactline,[lim(1) lim(2)],'-.'); hold on

figure(7)
plot(t+tplot,Menergy);title('Mechanical Energy');hold on;lim=ylim;plot(impactline,[lim(1) lim(2)],'-.'); hold on

energy_loss(i)=Menergy(size(Menergy,1))-fcn_energy(XX(:,i+1));
figure(8)
plot(i,energy_loss(i),'*');title('Energy Loss');hold on

step_length(i,:)=fcn_position_swingfoot(ye);
average_speed(i)=step_length(i,1)/te;

if i>=2
    accel(i)=(average_speed(i)-average_speed(i-1))/te;
end
tplot=tplot+te;



%thetadot_minus=c*ye(6:10)';
%checkpoint= betta(:,7)- H0*ye(6:10)'*(theta_minus-theta_plus)/thetadot_minus/M   
%betta

end
XX
timpact(i+1)=tplot;

figure (9)
%subplot(2,1,1);plot(timpact(2:end),average_speed,'rs');xlim([0 tplot]);xlabel('Time (s)','fontsize',16);ylabel('Average Speed (m/s)','fontsize',16);%title('Average Speed');
%subplot(3,1,2);plot(timpact(2:end),accel,'--rs');xlim([0 tplot]);title('Average Accelaration');hold on;% stairs(timpact(2:end),a_des,'-ok');hold on
%subplot(2,1,2);plot(tF,Fext,'--rs');ylim([0 12]);xlim([0 tplot]);xlabel('Time (s)','fontsize',16);ylabel('External Force','fontsize',16);%title('External Force');
plot(timpact(2:end),average_speed,'rs');xlim([0 tplot]);xlabel('Time (s)','fontsize',16);ylabel('Average Speed (m/s)','fontsize',16);%title('Average Speed');
figure (10)
plot(tF,Fext,'--rs');ylim([0 12]);xlim([0 tplot]);xlabel('Time (s)','fontsize',16);ylabel('External Force','fontsize',16);%title('External Force');

mean_accel=(average_speed(i)-average_speed(1))/(tplot-timpact(2))
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
