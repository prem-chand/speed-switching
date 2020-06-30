function cost=accel_cost(Fp)

M=6;
accel(1)=0;
te=0;
tplot=0;timpact=0;
a_des=[0 0.1 0.1 0.1 0.1 0.1 0.1];
clear  theta_minus theta_plus betta XX
global theta_minus theta_plus betta event

fix=open('/Users/mohamadshafiee/Desktop/Reasearch/programming/5 link RABBIT/fixedpointforfivelink/fixed_point_9.mat');
ival=fix.fpStates;

Fext=[0 0 0 Fp(1:5)];
%tF=[0 0.224 0.448 0.448+Fp(6) 0.448+Fp(6)+Fp(7) 0.448+Fp(6)+Fp(7)+Fp(8) 0.448+Fp(6)+Fp(7)+Fp(8)+Fp(9) 0.448+Fp(6)+Fp(7)+Fp(8)+Fp(9)+Fp(10)];
tF=[0,0.224051700858414,0.448114997580434,0.672186246208195,0.896265325583177,1.12035190722882,1.34444580807905,1.65];

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

end

cost=norm(a_des-accel);
