function cost=fixed_point_step(initial_alpha_red,Fext,tF)

M=6;lt=0.4;lf=0.4;
clear  theta_minus theta_plus betta XX
global theta_minus theta_plus betta event
%syms q5 real

q=initial_alpha_red(1:4);

%q5=solve(lt*cos(pi-q(1)-q(2)-q(4))+lf*cos(pi-q(1)-q(2))-lf*cos(q(1)+q(3)-pi)-lt*cos(q(1)+q(3)+q5-pi)+0.1)
ang=(lt*cos(pi-q(1)-q(2)-q(4))+lf*cos(pi-q(1)-q(2))-lf*cos(q(1)+q(3)-pi)+0.1)/lt
if (abs(ang)>1) 
    cost=10000;
else
q5=acos(ang)-q(1)-q(3)+pi

initial_alpha_red=[q q5(1) initial_alpha_red(5:end)];   % for GA
XX=initial_alpha_red(1:10)';

[betta,theta_minus,theta_plus]=fcn_alpha_red(initial_alpha_red');
event=XX(2);

options = odeset('Events',@touchdown5,'RelTol',1e-4,'AbsTol',1e-3);
[t,y,te,ye,ie] = ode45(@(t,Z) return_map5(t,Z,Fext,tF),[0 20],XX,options);



%if te==[]
%if y(size(y,1),3)-XX(2,1)~=0
if  size(te,1)==0
    cost=20000 ;
else
    %step_length=fcn_position_swingfoot(ye);
    %if step_length(1)<0 
    %    cost=20000
    %else
    %average_speed=step_length(1)/te
    xplus=impact_map5(ye);
    cost=norm(xplus-XX);%*.95+(1.1-average_speed)*0.05
end
% xplus=impact_map5(ye);
% cost=norm(xplus-XX)
%xplus=[ye(1);ye(3);ye(2);ye(5);ye(4)];
%cost=norm(xplus-XX(1:5));

end
end