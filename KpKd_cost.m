function cost=KpKd_cost(k)

M=6;
accel(1)=0;
te=0;
tplot=0;timpact=0;
F_step=0;
a_des=[0 0.1 0.1 0.1 0.1 0.1 0.1];
clear  theta_minus theta_plus betta XX
global theta_minus theta_plus betta event

fix=open('/Users/mohamadshafiee/Desktop/Reasearch/programming/5 link RABBIT/fixedpointforfivelink/fixed_point_9.mat');
ival=fix.fpStates;

ctr_index=15;
tF=[0 50];
Fext=[0  0];

[betta,theta_minus,theta_plus]=fcn_alpha_red(ival);
event=ival(2);
XX=ival(1:10);

for i=1:7

options = odeset('Events',@touchdown5,'RelTol',1e-5,'AbsTol',1e-5);
[t,y,te,ye,ie] = ode45(@(t,Z) return_map5(t+tplot,Z,Fext,tF),[0 20],XX(:,i),options);
XX(:,i+1)=impact_map5(ye);

step_length(i,:)=fcn_position_swingfoot(ye);
average_speed(i)=step_length(i,1)/te;


if i>=2
    accel(i)=(average_speed(i)-average_speed(i-1))/te;
end
tplot=tplot+te;

if i==1 
    F_step=k(1)*(a_des(i)-accel(i))+F_step;
elseif (i>=ctr_index+3) | (i<ctr_index)
    F_step=k(1)*(a_des(i)-accel(i))+k(2)*(accel(i)-accel(i-1))+F_step;
end

if a_des(i)==0 
    F_step=0;
end

% step-like force
Fext(2*i+1)=F_step;Fext(2*i+2)=F_step;
tF(2*i)=tplot;tF(2*i+1)=tplot+0.0001;tF(2*i+2)=tplot+3*te;

end

cost=norm(a_des-accel);
