clear;
close all;
clc;
M=6;
te=0;
tplot=0;timpact=0;
ntouch=0;ptouch=0;
accel(1)=0;
simout.t=[];simout.stk=[];simout.h=[]; simout.swk=[]; simout.swf=[]; simout.head=[];
tF=[0 1000];
% Fext=[1 1];    
Fext=[10 10];    
% fix=open('fixedpointforfivelink/fixed_point_9_gz_duanyi.mat');
fix=open('fixedpointforfivelink/fixed_point_vert2.mat');
% ival=fix.fpStates';
ival=fix.x;

clear  theta_minus theta_plus betta XX event
global theta_minus theta_plus betta event M_bez
M_bez=6;
XX=[];
XX=ival(1:10)';
event=XX(2);

c=[1 1 0 0.5 0];
[betta,theta_minus,theta_plus]=fcn_alpha_red(ival');

for i=1:5
    clear s
    theta_minus=XX(1,i)+XX(3,i)+0.5*XX(5,i); %c*q-
    theta_plus=c*XX(1:5,i);                  %c*q+
    deltheta=theta_plus - theta_minus;
    options = odeset('Events',@touchdown5,'RelTol',1e-5,'AbsTol',1e-4);
    [t,y,te,ye,ie] = ode15s(@(t,Z) return_map5(t+tplot,Z,Fext,tF,betta,theta_minus,theta_plus,M_bez),[0:0.001:20],XX(:,i),options);
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

    [x_plus,Fground]=impact_map5(ye);
    XX(:,i+1)=x_plus;
    gam0_plus=gamma0(XX(:,i+1));
    gam0_minus=gamma0(ye);
    delta_zero=(gam0_plus*XX(6:10,i+1))/(gam0_minus*ye(6:10)');
    timpact(i)=tplot;

    impactline=[timpact(i) timpact(i)];

    clear u h hdot GRF s k1 k2 ksi1 ksi2 Vzero Vzero_ksi1 theta Wf Fvir;

    for tp=1:size(t)
        [fx,gx,gforce]=fcn_stance_dynamics5(y(tp,:));
        Fext_plot = [interp1(tF,Fext,t(tp)+tplot);0];
        theta(tp)=c*y(tp,1:5)';
        s(tp)=(theta(tp)-theta_plus)/(theta_minus-theta_plus);
        k1(tp)=c*([fcn_Dh_Dq(s(tp),betta,deltheta) ; gamma0(y(tp,1:5))]\[zeros(4,1) ; 1]);
        k2(tp)=kappa2(y(tp,1:5));
        k3=kappa3(y(tp,1:5));
        Fvir(tp)=(k3*Fext_plot)/k1(tp);
        ksi1(tp)=theta(tp);
        ksi2(tp)=gamma0(y(tp,1:5)')*y(tp,6:10)';
%         if M_bez==6
            [h(tp,:),hdot(tp,:)]=fcn_liederivative(y(tp,:),betta,theta_minus,theta_plus,M_bez);
%         elseif M_bez==5
%             [h(tp,:),hdot(tp,:)]=fcn_liederivative_m5(y(tp,:),betta);
%         end
    end
    Vzero=trapz(theta,-k2./k1-Fvir);
    Wf=trapz(theta,Fvir);
    
    for j=2:1:length(t)
        Vzero_ksi1(j)=trapz(theta(1:j),-k2(1:j)./k1(1:j)-Fvir(1:j));
    end
    K=max(Vzero_ksi1);
    
    ksi2_new=sqrt(2*(0.5*19.65^2+Wf/(1-delta_zero^2)));
    
%     figure(1)
%     subplot(1,2,1); plot(theta,k1); hold on;
%     subplot(1,2,2); plot(theta,k2); hold on;
    Vzero_ksi1;
    Vzero;
    boa=sqrt(2*K/delta_zero^2);
    % ksi2(end)
    figure(2)
%     plot(i,Vzero,'o'); hold on;
%     plot(i,ksi2(end),'o'); hold on;
%     plot(theta,c*y(:,6:10)'); hold on;
%     plot(ksi1,ksi2); hold on;
%     plot([theta_plus theta_minus],[boa boa]); hold on;
%     (delta_zero^2/(1-delta_zero^2))*Vzero+K

    step_length(i,:)=fcn_position_swingfoot(ye);
    average_speed(i)=step_length(i,1)/te;
    
    if i>0
    output_error=[h hdot];
    norm_output_error = (sqrt(sum(output_error.^2,2)));
    plot(t+tplot,norm_output_error); hold on;
    end

    if i>=2
        accel(i)=(average_speed(i)-average_speed(i-1))/te;
    end
    tplot=tplot+te;

end
mean(accel)
timpact(i+1)=tplot;
mean_accel=(average_speed(i)-average_speed(1))/(tplot-timpact(2));
hold off;