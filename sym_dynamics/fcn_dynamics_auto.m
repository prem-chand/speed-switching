function [G,D,C,J,v2v] = fcn_dynamics_auto(q)

  model_params;

  G(1,1)=g*(Mt*(sin(q(1) + q(2) + q(4))*(lt - pMt) + lf*sin(q(1) + q(2)) - lf*sin(q(1) + q(3)) + lt*...
         sin(q(1) + q(2) + q(4)) - pMt*sin(q(1) + q(3) + q(5))) + MT*(lf*sin(q(1) + q(2)) - pMT*sin(q(1)) + lt*...
         sin(q(1) + q(2) + q(4))) + Mf*(lf*sin(q(1) + q(2)) - pMf*sin(q(1) + q(3)) + sin(q(1) + q(2))*(lf - pMf) + 2*...
         lt*sin(q(1) + q(2) + q(4))));
  G(2,1)=g*(Mf*(lf*sin(q(1) + q(2)) + sin(q(1) + q(2))*(lf - pMf) + 2*lt*sin(q(1) + q(2) + q(4))) +...
          MT*(lf*sin(q(1) + q(2)) + lt*sin(q(1) + q(2) + q(4))) + Mt*(sin(q(1) + q(2) + q(4))*(lt - pMt) + lf*...
         sin(q(1) + q(2)) + lt*sin(q(1) + q(2) + q(4))));
  G(3,1)=-g*(Mt*(lf*sin(q(1) + q(3)) + pMt*sin(q(1) + q(3) + q(5))) + Mf*pMf*sin(q(1) + q(3)));
  G(4,1)=g*sin(q(1) + q(2) + q(4))*(MT*lt + 2*Mf*lt + 2*Mt*lt - Mt*pMt);
  G(5,1)=-Mt*g*pMt*sin(q(1) + q(3) + q(5));

   D(1,1)=IT + 2*If + 2*It + MT*lf^2 + MT*lt^2 + 2*Mf*lf^2 + 2*Mf*lt^2 + 2*Mt*lf^2 + 2*Mt*lt^2 + MT*...
         pMT^2 + 2*Mf*pMf^2 + 2*Mt*pMt^2 - 2*Mt*lf^2*cos(q(2) - q(3)) - 2*Mf*lf*pMf - 2*Mt*lt*pMt - 2*MT*lt*pMT*...
         cos(q(2) + q(4)) + 2*MT*lf*lt*cos(q(4)) + 4*Mf*lf*lt*cos(q(4)) + 2*Mt*lf*lt*cos(q(4)) - 2*MT*lf*pMT*...
         cos(q(2)) - 2*Mf*lt*pMf*cos(q(4)) + 2*Mt*lf*pMt*cos(q(5)) - 2*Mt*lf*pMt*cos(q(2) - q(3) - q(5)) - 2*Mf*lf*pMf*...
         cos(q(2) - q(3)) - 2*Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5)) - 2*Mt*lf*lt*cos(q(2) - q(3) + q(4)) - 2*Mf*lt*pMf*...
         cos(q(2) - q(3) + q(4));
  D(1,2)=If + It + MT*lf^2 + MT*lt^2 + 2*Mf*lf^2 + 2*Mf*lt^2 + Mt*lf^2 + 2*Mt*lt^2 + Mf*pMf^2 + Mt*...
         pMt^2 - Mt*lf^2*cos(q(2) - q(3)) - 2*Mf*lf*pMf - 2*Mt*lt*pMt - MT*lt*pMT*cos(q(2) + q(4)) + 2*MT*lf*lt*...
         cos(q(4)) + 4*Mf*lf*lt*cos(q(4)) + 2*Mt*lf*lt*cos(q(4)) - MT*lf*pMT*cos(q(2)) - 2*Mf*lt*pMf*cos(q(4)) - Mt*lf*...
         pMt*cos(q(2) - q(3) - q(5)) - Mf*lf*pMf*cos(q(2) - q(3)) - Mt*lt*pMt*cos(q(2) - q(3) +...
          q(4) - q(5)) - Mt*lf*lt*cos(q(2) - q(3) + q(4)) - Mf*lt*pMf*cos(q(2) - q(3) + q(4));
  D(1,3)=If + It + Mt*lf^2 + Mf*pMf^2 + Mt*pMt^2 - Mt*lf^2*cos(q(2) - q(3)) + 2*Mt*lf*pMt*...
         cos(q(5)) - Mt*lf*pMt*cos(q(2) - q(3) - q(5)) - Mf*lf*pMf*cos(q(2) - q(3)) - Mt*lt*pMt*cos(q(2) - q(3) +...
          q(4) - q(5)) - Mt*lf*lt*cos(q(2) - q(3) + q(4)) - Mf*lt*pMf*cos(q(2) - q(3) + q(4));
  D(1,4)=It + MT*lt^2 + 2*Mf*lt^2 + 2*Mt*lt^2 + Mt*pMt^2 - 2*Mt*lt*pMt - MT*lt*pMT*cos(q(2) +...
          q(4)) + MT*lf*lt*cos(q(4)) + 2*Mf*lf*lt*cos(q(4)) + Mt*lf*lt*cos(q(4)) - Mf*lt*pMf*cos(q(4)) - Mt*lt*pMt*...
         cos(q(2) - q(3) + q(4) - q(5)) - Mt*lf*lt*cos(q(2) - q(3) + q(4)) - Mf*lt*pMf*cos(q(2) - q(3) + q(4));
  D(1,5)=It + Mt*pMt^2 + Mt*lf*pMt*cos(q(5)) - Mt*lf*pMt*cos(q(2) - q(3) - q(5)) - Mt*lt*pMt*...
         cos(q(2) - q(3) + q(4) - q(5));
  D(2,1)=If + It + MT*lf^2 + MT*lt^2 + 2*Mf*lf^2 + 2*Mf*lt^2 + Mt*lf^2 + 2*Mt*lt^2 + Mf*pMf^2 + Mt*...
         pMt^2 - Mt*lf^2*cos(q(2) - q(3)) - 2*Mf*lf*pMf - 2*Mt*lt*pMt - MT*lt*pMT*cos(q(2) + q(4)) + 2*MT*lf*lt*...
         cos(q(4)) + 4*Mf*lf*lt*cos(q(4)) + 2*Mt*lf*lt*cos(q(4)) - MT*lf*pMT*cos(q(2)) - 2*Mf*lt*pMf*cos(q(4)) - Mt*lf*...
         pMt*cos(q(2) - q(3) - q(5)) - Mf*lf*pMf*cos(q(2) - q(3)) - Mt*lt*pMt*cos(q(2) - q(3) +...
          q(4) - q(5)) - Mt*lf*lt*cos(q(2) - q(3) + q(4)) - Mf*lt*pMf*cos(q(2) - q(3) + q(4));
  D(2,2)=If + It + MT*lf^2 + MT*lt^2 + 2*Mf*lf^2 + 2*Mf*lt^2 + Mt*lf^2 + 2*Mt*lt^2 + Mf*pMf^2 + Mt*...
         pMt^2 - 2*Mf*lf*pMf - 2*Mt*lt*pMt + 2*MT*lf*lt*cos(q(4)) + 4*Mf*lf*lt*cos(q(4)) + 2*Mt*lf*lt*cos(q(4)) - 2*Mf*...
         lt*pMf*cos(q(4));
  D(2,3)=- Mt*lf^2*cos(q(2) - q(3)) - Mt*lf*pMt*cos(q(2) - q(3) - q(5)) - Mf*lf*pMf*...
         cos(q(2) - q(3)) - Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5)) - Mt*lf*lt*cos(q(2) - q(3) + q(4)) - Mf*lt*pMf*...
         cos(q(2) - q(3) + q(4));
  D(2,4)=It + MT*lt^2 + 2*Mf*lt^2 + 2*Mt*lt^2 + Mt*pMt^2 - 2*Mt*lt*pMt + MT*lf*lt*cos(q(4)) + 2*Mf*...
         lf*lt*cos(q(4)) + Mt*lf*lt*cos(q(4)) - Mf*lt*pMf*cos(q(4));
  D(2,5)=- Mt*lf*pMt*cos(q(2) - q(3) - q(5)) - Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5));
  D(3,1)=If + It + Mt*lf^2 + Mf*pMf^2 + Mt*pMt^2 - Mt*lf^2*cos(q(2) - q(3)) + 2*Mt*lf*pMt*...
         cos(q(5)) - Mt*lf*pMt*cos(q(2) - q(3) - q(5)) - Mf*lf*pMf*cos(q(2) - q(3)) - Mt*lt*pMt*cos(q(2) - q(3) +...
          q(4) - q(5)) - Mt*lf*lt*cos(q(2) - q(3) + q(4)) - Mf*lt*pMf*cos(q(2) - q(3) + q(4));
  D(3,2)=- Mt*lf^2*cos(q(2) - q(3)) - Mt*lf*pMt*cos(q(2) - q(3) - q(5)) - Mf*lf*pMf*...
         cos(q(2) - q(3)) - Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5)) - Mt*lf*lt*cos(q(2) - q(3) + q(4)) - Mf*lt*pMf*...
         cos(q(2) - q(3) + q(4));
  D(3,3)=If + It + Mt*lf^2 + Mf*pMf^2 + Mt*pMt^2 + 2*Mt*lf*pMt*cos(q(5));
  D(3,4)=- Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5)) - Mt*lf*lt*cos(q(2) - q(3) + q(4)) - Mf*lt*pMf*...
         cos(q(2) - q(3) + q(4));
  D(3,5)=It + Mt*pMt^2 + Mt*lf*pMt*cos(q(5));
  D(4,1)=It + MT*lt^2 + 2*Mf*lt^2 + 2*Mt*lt^2 + Mt*pMt^2 - 2*Mt*lt*pMt - MT*lt*pMT*cos(q(2) +...
          q(4)) + MT*lf*lt*cos(q(4)) + 2*Mf*lf*lt*cos(q(4)) + Mt*lf*lt*cos(q(4)) - Mf*lt*pMf*cos(q(4)) - Mt*lt*pMt*...
         cos(q(2) - q(3) + q(4) - q(5)) - Mt*lf*lt*cos(q(2) - q(3) + q(4)) - Mf*lt*pMf*cos(q(2) - q(3) + q(4));
  D(4,2)=It + MT*lt^2 + 2*Mf*lt^2 + 2*Mt*lt^2 + Mt*pMt^2 - 2*Mt*lt*pMt + MT*lf*lt*cos(q(4)) + 2*Mf*...
         lf*lt*cos(q(4)) + Mt*lf*lt*cos(q(4)) - Mf*lt*pMf*cos(q(4));
  D(4,3)=- Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5)) - Mt*lf*lt*cos(q(2) - q(3) + q(4)) - Mf*lt*pMf*...
         cos(q(2) - q(3) + q(4));
  D(4,4)=It + MT*lt^2 + 2*Mf*lt^2 + Mt*lt^2 + Mt*(lt - pMt)^2;
  D(4,5)=-Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5));
  D(5,1)=It + Mt*pMt^2 + Mt*lf*pMt*cos(q(5)) - Mt*lf*pMt*cos(q(2) - q(3) - q(5)) - Mt*lt*pMt*...
         cos(q(2) - q(3) + q(4) - q(5));
  D(5,2)=- Mt*lf*pMt*cos(q(2) - q(3) - q(5)) - Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5));
  D(5,3)=It + Mt*pMt^2 + Mt*lf*pMt*cos(q(5));
  D(5,4)=-Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5));
  D(5,5)=It + Mt*pMt^2;

   C(1,1)=Mt*dq(2)*lf^2*sin(q(2) - q(3)) - Mt*dq(3)*lf^2*sin(q(2) - q(3)) + Mf*dq(2)*lf*pMf*...
         sin(q(2) - q(3)) - Mf*dq(3)*lf*pMf*sin(q(2) - q(3)) + Mt*dq(2)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5)) - Mt*dq(3)*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5)) + Mt*dq(4)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5)) - Mt*dq(5)*lt*pMt*sin(q(2) - q(3) +...
          q(4) - q(5)) + Mt*dq(2)*lf*lt*sin(q(2) - q(3) + q(4)) - Mt*dq(3)*lf*lt*sin(q(2) - q(3) + q(4)) + Mt*dq(4)*lf*lt*...
         sin(q(2) - q(3) + q(4)) + Mf*dq(2)*lt*pMf*sin(q(2) - q(3) + q(4)) - Mf*dq(3)*lt*pMf*sin(q(2) - q(3) + q(4)) + Mf*...
         dq(4)*lt*pMf*sin(q(2) - q(3) + q(4)) + MT*dq(2)*lt*pMT*sin(q(2) + q(4)) + MT*dq(4)*lt*pMT*sin(q(2) +...
          q(4)) - MT*dq(4)*lf*lt*sin(q(4)) - 2*Mf*dq(4)*lf*lt*sin(q(4)) - Mt*dq(4)*lf*lt*sin(q(4)) + MT*dq(2)*lf*pMT*...
         sin(q(2)) + Mf*dq(4)*lt*pMf*sin(q(4)) - Mt*dq(5)*lf*pMt*sin(q(5)) + Mt*dq(2)*lf*pMt*...
         sin(q(2) - q(3) - q(5)) - Mt*dq(3)*lf*pMt*sin(q(2) - q(3) - q(5)) - Mt*dq(5)*lf*pMt*sin(q(2) - q(3) - q(5));
  C(1,2)=dq(1)*(Mt*lf^2*sin(q(2) - q(3)) + Mt*lf*lt*sin(q(2) - q(3) + q(4)) + Mf*lt*pMf*...
         sin(q(2) - q(3) + q(4)) + MT*lt*pMT*sin(q(2) + q(4)) + MT*lf*pMT*sin(q(2)) + Mt*lf*pMt*sin(q(2) - q(3) - q(5)) + Mf*...
         lf*pMf*sin(q(2) - q(3)) + Mt*lt*pMt*sin(q(2) - q(3) + q(4) - q(5))) + dq(2)*(Mt*lf^2*sin(q(2) - q(3)) +...
          Mt*lf*lt*sin(q(2) - q(3) + q(4)) + Mf*lt*pMf*sin(q(2) - q(3) + q(4)) + MT*lt*pMT*sin(q(2) + q(4)) + MT*...
         lf*pMT*sin(q(2)) + Mt*lf*pMt*sin(q(2) - q(3) - q(5)) + Mf*lf*pMf*sin(q(2) - q(3)) + Mt*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5))) + dq(4)*lt*(Mf*pMf*sin(q(4)) + Mt*pMt*sin(q(2) - q(3) + q(4) - q(5)) + Mt*lf*...
         sin(q(2) - q(3) + q(4)) + Mf*pMf*sin(q(2) - q(3) + q(4)) + MT*pMT*sin(q(2) + q(4)) - MT*lf*sin(q(4)) - 2*Mf*lf*...
         sin(q(4)) - Mt*lf*sin(q(4)));
  C(1,3)=- dq(1)*(Mt*lf^2*sin(q(2) - q(3)) + Mt*lf*lt*sin(q(2) - q(3) + q(4)) + Mf*lt*pMf*...
         sin(q(2) - q(3) + q(4)) + Mt*lf*pMt*sin(q(2) - q(3) - q(5)) + Mf*lf*pMf*sin(q(2) - q(3)) + Mt*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5))) - dq(3)*(Mt*lf^2*sin(q(2) - q(3)) + Mt*lf*lt*sin(q(2) - q(3) + q(4)) + Mf*lt*pMf*...
         sin(q(2) - q(3) + q(4)) + Mt*lf*pMt*sin(q(2) - q(3) - q(5)) + Mf*lf*pMf*sin(q(2) - q(3)) + Mt*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5))) - Mt*dq(5)*pMt*(lt*sin(q(2) - q(3) + q(4) - q(5)) + lf*sin(q(5)) + lf*sin(q(2) - q(3) - q(5)));
  C(1,4)=lt*(dq(1) + dq(2) + dq(4))*(Mf*pMf*sin(q(4)) + Mt*pMt*sin(q(2) - q(3) + q(4) - q(5)) + Mt*...
         lf*sin(q(2) - q(3) + q(4)) + Mf*pMf*sin(q(2) - q(3) + q(4)) + MT*pMT*sin(q(2) + q(4)) - MT*lf*...
         sin(q(4)) - 2*Mf*lf*sin(q(4)) - Mt*lf*sin(q(4)));
  C(1,5)=-Mt*pMt*(dq(1) + dq(3) + dq(5))*(lt*sin(q(2) - q(3) + q(4) - q(5)) + lf*sin(q(5)) + lf*sin(q(2) - q(3) - q(5)));
  C(2,1)=Mf*dq(4)*lt*pMf*sin(q(4)) - Mt*dq(3)*lf^2*sin(q(2) - q(3)) - Mf*dq(1)*lf*pMf*...
         sin(q(2) - q(3)) - Mf*dq(3)*lf*pMf*sin(q(2) - q(3)) - Mt*dq(1)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5)) - Mt*dq(3)*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5)) - Mt*dq(5)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5)) - Mt*dq(1)*lf*lt*sin(q(2) - q(3) +...
          q(4)) - Mt*dq(3)*lf*lt*sin(q(2) - q(3) + q(4)) - Mf*dq(1)*lt*pMf*sin(q(2) - q(3) + q(4)) - Mf*dq(3)*lt*pMf*...
         sin(q(2) - q(3) + q(4)) - MT*dq(1)*lt*pMT*sin(q(2) + q(4)) - MT*dq(4)*lf*lt*sin(q(4)) - 2*Mf*dq(4)*lf*lt*...
         sin(q(4)) - Mt*dq(4)*lf*lt*sin(q(4)) - MT*dq(1)*lf*pMT*sin(q(2)) - Mt*dq(1)*lf^2*sin(q(2) - q(3)) - Mt*dq(1)*lf*...
         pMt*sin(q(2) - q(3) - q(5)) - Mt*dq(3)*lf*pMt*sin(q(2) - q(3) - q(5)) - Mt*dq(5)*lf*pMt*sin(q(2) - q(3) - q(5));
  C(2,2)=-dq(4)*lt*sin(q(4))*(MT*lf + 2*Mf*lf + Mt*lf - Mf*pMf);
  C(2,3)=- dq(1)*(Mt*lf^2*sin(q(2) - q(3)) + Mt*lf*lt*sin(q(2) - q(3) + q(4)) + Mf*lt*pMf*...
         sin(q(2) - q(3) + q(4)) + Mt*lf*pMt*sin(q(2) - q(3) - q(5)) + Mf*lf*pMf*sin(q(2) - q(3)) + Mt*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5))) - dq(3)*(Mt*lf^2*sin(q(2) - q(3)) + Mt*lf*lt*sin(q(2) - q(3) + q(4)) + Mf*lt*pMf*...
         sin(q(2) - q(3) + q(4)) + Mt*lf*pMt*sin(q(2) - q(3) - q(5)) + Mf*lf*pMf*sin(q(2) - q(3)) + Mt*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5))) - Mt*dq(5)*pMt*(lt*sin(q(2) - q(3) + q(4) - q(5)) + lf*sin(q(2) - q(3) - q(5)));
  C(2,4)=-lt*sin(q(4))*(dq(1) + dq(2) + dq(4))*(MT*lf + 2*Mf*lf + Mt*lf - Mf*pMf);
  C(2,5)=-Mt*pMt*(lt*sin(q(2) - q(3) + q(4) - q(5)) + lf*sin(q(2) - q(3) - q(5)))*(dq(1) + dq(3) + dq(5));
  C(3,1)=Mt*dq(1)*lf^2*sin(q(2) - q(3)) + Mt*dq(2)*lf^2*sin(q(2) - q(3)) + Mf*dq(1)*lf*pMf*...
         sin(q(2) - q(3)) + Mf*dq(2)*lf*pMf*sin(q(2) - q(3)) + Mt*dq(1)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5)) + Mt*dq(2)*lt*...
         pMt*sin(q(2) - q(3) + q(4) - q(5)) + Mt*dq(4)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5)) + Mt*dq(1)*lf*lt*...
         sin(q(2) - q(3) + q(4)) + Mt*dq(2)*lf*lt*sin(q(2) - q(3) + q(4)) + Mt*dq(4)*lf*lt*sin(q(2) - q(3) + q(4)) + Mf*dq(1)*...
         lt*pMf*sin(q(2) - q(3) + q(4)) + Mf*dq(2)*lt*pMf*sin(q(2) - q(3) + q(4)) + Mf*dq(4)*lt*pMf*...
         sin(q(2) - q(3) + q(4)) - Mt*dq(5)*lf*pMt*sin(q(5)) + Mt*dq(1)*lf*pMt*sin(q(2) - q(3) - q(5)) + Mt*dq(2)*lf*pMt*sin(q(2) - q(3) - q(5));
  C(3,2)=dq(1)*(Mt*lf^2*sin(q(2) - q(3)) + Mt*lf*lt*sin(q(2) - q(3) + q(4)) + Mf*lt*pMf*...
         sin(q(2) - q(3) + q(4)) + Mt*lf*pMt*sin(q(2) - q(3) - q(5)) + Mf*lf*pMf*sin(q(2) - q(3)) + Mt*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5))) + dq(2)*(Mt*lf^2*sin(q(2) - q(3)) + Mt*lf*lt*sin(q(2) - q(3) + q(4)) + Mf*lt*pMf*...
         sin(q(2) - q(3) + q(4)) + Mt*lf*pMt*sin(q(2) - q(3) - q(5)) + Mf*lf*pMf*sin(q(2) - q(3)) + Mt*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5))) + dq(4)*lt*(Mt*pMt*sin(q(2) - q(3) + q(4) - q(5)) + Mt*lf*sin(q(2) - q(3) + q(4)) +...
          Mf*pMf*sin(q(2) - q(3) + q(4)));
  C(3,3)=-Mt*dq(5)*lf*pMt*sin(q(5));
  C(3,4)=lt*(Mt*pMt*sin(q(2) - q(3) + q(4) - q(5)) + Mt*lf*sin(q(2) - q(3) + q(4)) + Mf*pMf*...
         sin(q(2) - q(3) + q(4)))*(dq(1) + dq(2) + dq(4));
  C(3,5)=-Mt*lf*pMt*sin(q(5))*(dq(1) + dq(3) + dq(5));
  C(4,1)=MT*dq(1)*lf*lt*sin(q(4)) - Mt*dq(3)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5)) - Mt*dq(5)*lt*...
         pMt*sin(q(2) - q(3) + q(4) - q(5)) - Mt*dq(1)*lf*lt*sin(q(2) - q(3) + q(4)) - Mt*dq(3)*lf*lt*...
         sin(q(2) - q(3) + q(4)) - Mf*dq(1)*lt*pMf*sin(q(2) - q(3) + q(4)) - Mf*dq(3)*lt*pMf*sin(q(2) - q(3) + q(4)) - MT*...
         dq(1)*lt*pMT*sin(q(2) + q(4)) - Mt*dq(1)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5)) + MT*dq(2)*lf*lt*...
         sin(q(4)) + 2*Mf*dq(1)*lf*lt*sin(q(4)) + 2*Mf*dq(2)*lf*lt*sin(q(4)) + Mt*dq(1)*lf*lt*sin(q(4)) + Mt*dq(2)*lf*...
         lt*sin(q(4)) - Mf*dq(1)*lt*pMf*sin(q(4)) - Mf*dq(2)*lt*pMf*sin(q(4));
  C(4,2)=lt*sin(q(4))*(dq(1) + dq(2))*(MT*lf + 2*Mf*lf + Mt*lf - Mf*pMf);
  C(4,3)=- dq(1)*lt*(Mt*pMt*sin(q(2) - q(3) + q(4) - q(5)) + Mt*lf*sin(q(2) - q(3) + q(4)) + Mf*pMf*...
         sin(q(2) - q(3) + q(4))) - dq(3)*lt*(Mt*pMt*sin(q(2) - q(3) + q(4) - q(5)) + Mt*lf*sin(q(2) - q(3) + q(4)) + Mf*pMf*...
         sin(q(2) - q(3) + q(4))) - Mt*dq(5)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5));
  C(4,4)=0;
  C(4,5)=-Mt*lt*pMt*sin(q(2) - q(3) + q(4) - q(5))*(dq(1) + dq(3) + dq(5));
  C(5,1)=Mt*pMt*(dq(1)*lf*sin(q(5)) + dq(3)*lf*sin(q(5)) + dq(1)*lf*sin(q(2) - q(3) - q(5)) + dq(2)*...
         lf*sin(q(2) - q(3) - q(5)) + dq(1)*lt*sin(q(2) - q(3) + q(4) - q(5)) + dq(2)*lt*sin(q(2) - q(3) +...
          q(4) - q(5)) + dq(4)*lt*sin(q(2) - q(3) + q(4) - q(5)));
  C(5,2)=dq(1)*(Mt*lf*pMt*sin(q(2) - q(3) - q(5)) + Mt*lt*pMt*sin(q(2) - q(3) + q(4) - q(5))) +...
          dq(2)*(Mt*lf*pMt*sin(q(2) - q(3) - q(5)) + Mt*lt*pMt*sin(q(2) - q(3) + q(4) - q(5))) + Mt*dq(4)*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5));
  C(5,3)=Mt*lf*pMt*sin(q(5))*(dq(1) + dq(3));
  C(5,4)=Mt*lt*pMt*sin(q(2) - q(3) + q(4) - q(5))*(dq(1) + dq(2) + dq(4));
  C(5,5)=0;

   J(1,1)=lT*cos(q(1)) - lf*cos(q(1) + q(2)) - lt*cos(q(1) + q(2) + q(4));
  J(1,2)=- lf*cos(q(1) + q(2)) - lt*cos(q(1) + q(2) + q(4));
  J(1,3)=0;
  J(1,4)=-lt*cos(q(1) + q(2) + q(4));
  J(1,5)=0;
  J(2,1)=lf*sin(q(1) + q(2)) - lT*sin(q(1)) + lt*sin(q(1) + q(2) + q(4));
  J(2,2)=lf*sin(q(1) + q(2)) + lt*sin(q(1) + q(2) + q(4));
  J(2,3)=0;
  J(2,4)=lt*sin(q(1) + q(2) + q(4));
  J(2,5)=0;

   v2v(1,1)=dq(3)*(lf*cos(q(1) + q(3)) + lt*cos(q(1) + q(3) + q(5))) - dq(2)*(lf*cos(q(1) + q(2)) +...
          lt*cos(q(1) + q(2) + q(4))) - dq(1)*(lf*cos(q(1) + q(2)) - lf*cos(q(1) + q(3)) + lt*cos(q(1) + q(2) +...
          q(4)) - lt*cos(q(1) + q(3) + q(5))) - dq(4)*lt*cos(q(1) + q(2) + q(4)) + dq(5)*lt*cos(q(1) + q(3) + q(5));
  v2v(2,1)=dq(1)*(lf*sin(q(1) + q(2)) - lf*sin(q(1) + q(3)) + lt*sin(q(1) + q(2) + q(4)) - lt*...
         sin(q(1) + q(3) + q(5))) + dq(2)*(lf*sin(q(1) + q(2)) + lt*sin(q(1) + q(2) + q(4))) - dq(3)*(lf*sin(q(1) +...
          q(3)) + lt*sin(q(1) + q(3) + q(5))) + dq(4)*lt*sin(q(1) + q(2) + q(4)) - dq(5)*lt*sin(q(1) + q(3) + q(5));

 