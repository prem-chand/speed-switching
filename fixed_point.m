function cost=fixed_point(alpha_red,Fext,tF)
global betta XX
betta(1,3:6)=alpha_red(1:4);
betta(2,3:6)=alpha_red(5:8);
betta(3,3:6)=alpha_red(9:12);
betta(4,3:6)=alpha_red(13:16);

options = odeset('Events',@touchdown5,'RelTol',1e-9,'AbsTol',1e-9);
[t,y,te,ye,ie] = ode45(@(t,Z) return_map5(t,Z,Fext,tF),[0 20],XX,options);
ye
te
xplus=impact_map5(ye);
cost=norm(xplus-XX)
%xplus=[ye(1);ye(3);ye(2);ye(5);ye(4)];
%cost=norm(xplus-XX(1:5));

end