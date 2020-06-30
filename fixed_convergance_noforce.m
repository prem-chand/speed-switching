function cost=fixed_convergance_noforce(initial_alpha_red,Fext,tF)

M=6;
clear  theta_minus theta_plus betta XX
global theta_minus theta_plus betta event XX_noforce

XX=XX_noforce(:,1) ;   % for GA
%XX=initial_alpha_red(1:10);
aug_initial_alpha=[XX' initial_alpha_red];

[betta,theta_minus,theta_plus]=fcn_alpha_red(aug_initial_alpha');
event=XX(2);

options = odeset('Events',@touchdown5,'RelTol',1e-5,'AbsTol',1e-5);
[t,y,te,ye,ie] = ode45(@(t,Z) return_map5(t,Z,Fext,tF),[0 20],XX,options);



%if te==[]
%if y(size(y,1),3)-XX(2,1)~=0
if  size(te,1)==0
    cost=20000 ;
else
 
 xplus=impact_map5(ye);
 cost=norm(xplus-XX_noforce(:,2));
%xplus=[ye(1);ye(3);ye(2);ye(5);ye(4)];
%cost=norm(xplus-XX(1:5));

end