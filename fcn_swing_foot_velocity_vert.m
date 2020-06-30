function v2v = fcn_swing_foot_velocity_vert(x)

q=x(1:5);
dq=x(6:10);

g=9.8;lT=0.63;lf=0.4;lt=0.4;MT=12;Mf=6.8;Mt=3.2;IT=1.33;If=0.47;It=0.2;pMT=0.24;pMf=0.11;pMt=0.24;

v2v(1,1)=dq(3)*(lf*cos(q(1) + q(3)) + lt*cos(q(1) + q(3) + q(5))) - dq(2)*(lf*cos(q(1) + q(2)) +...
          lt*cos(q(1) + q(2) + q(4))) - dq(1)*(lf*cos(q(1) + q(2)) - lf*cos(q(1) + q(3)) + lt*cos(q(1) + q(2) +...
          q(4)) - lt*cos(q(1) + q(3) + q(5))) - dq(4)*lt*cos(q(1) + q(2) + q(4)) + dq(5)*lt*cos(q(1) + q(3) + q(5));
  v2v(2,1)=dq(1)*(lf*sin(q(1) + q(2)) - lf*sin(q(1) + q(3)) + lt*sin(q(1) + q(2) + q(4)) - lt*...
         sin(q(1) + q(3) + q(5))) + dq(2)*(lf*sin(q(1) + q(2)) + lt*sin(q(1) + q(2) + q(4))) - dq(3)*(lf*sin(q(1) +...
          q(3)) + lt*sin(q(1) + q(3) + q(5))) + dq(4)*lt*sin(q(1) + q(2) + q(4)) - dq(5)*lt*sin(q(1) + q(3) + q(5));
end

