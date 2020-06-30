function [G,D,C] = fcn_dynamics_matrices2(q)

  model_params;

  G(1,1)=-conj(g)*(conj(Mt)*(sin(conj(q(1)) + conj(q(2)) + conj(q(4)))*conj(pMt) - 2*...
         sin(conj(q(1)) + conj(q(2)) + conj(q(4)))*conj(lt) + sin(conj(q(1)) + conj(q(3)) + conj(q(5)))*...
         conj(pMt) - sin(conj(q(1)) + conj(q(2)))*conj(lf) + sin(conj(q(1)) + conj(q(3)))*conj(lf)) - conj(MT)*(sin(conj(q(1)) +...
          conj(q(2)) + conj(q(4)))*conj(lt) - sin(conj(q(1)))*conj(pMT) + sin(conj(q(1)) + conj(q(2)))*conj(lf)) +...
          conj(Mfa)*(sin(conj(q(1)))*conj(psh) - sin(conj(q(1)) + conj(q(2)) + conj(q(4)))*conj(lt) - sin(conj(q(1)) +...
          conj(q(2)))*conj(lf) + sin(conj(q(1)) + conj(q(6)))*conj(lua) + sin(conj(q(1)) + conj(q(6)) - conj(q(7)))*...
         conj(pMfa)) + conj(Mua)*(sin(conj(q(1)))*conj(psh) - sin(conj(q(1)) + conj(q(2)) + conj(q(4)))*...
         conj(lt) - sin(conj(q(1)) + conj(q(2)))*conj(lf) + sin(conj(q(1)) + conj(q(6)))*conj(pMua)) - conj(Mf)*(2*sin(conj(q(1)) +...
          conj(q(2)) + conj(q(4)))*conj(lt) + 2*sin(conj(q(1)) + conj(q(2)))*conj(lf) - sin(conj(q(1)) + conj(q(2)))*...
         conj(pMf) - sin(conj(q(1)) + conj(q(3)))*conj(pMf)));
  G(2,1)=conj(g)*(conj(MT)*(sin(conj(q(1)) + conj(q(2)) + conj(q(4)))*conj(lt) + sin(conj(q(1)) +...
          conj(q(2)))*conj(lf)) + conj(Mfa)*(sin(conj(q(1)) + conj(q(2)) + conj(q(4)))*conj(lt) + sin(conj(q(1)) +...
          conj(q(2)))*conj(lf)) + conj(Mua)*(sin(conj(q(1)) + conj(q(2)) + conj(q(4)))*conj(lt) + sin(conj(q(1)) +...
          conj(q(2)))*conj(lf)) + conj(Mt)*(sin(conj(q(1)) + conj(q(2)) + conj(q(4)))*(conj(lt) - conj(pMt)) +...
          sin(conj(q(1)) + conj(q(2)) + conj(q(4)))*conj(lt) + sin(conj(q(1)) + conj(q(2)))*conj(lf)) + conj(Mf)*...
         (sin(conj(q(1)) + conj(q(2)))*(conj(lf) - conj(pMf)) + 2*sin(conj(q(1)) + conj(q(2)) + conj(q(4)))*conj(lt) +...
          sin(conj(q(1)) + conj(q(2)))*conj(lf)));
  G(3,1)=-conj(g)*(conj(Mt)*(sin(conj(q(1)) + conj(q(3)) + conj(q(5)))*conj(pMt) + sin(conj(q(1)) +...
          conj(q(3)))*conj(lf)) + sin(conj(q(1)) + conj(q(3)))*conj(Mf)*conj(pMf));
  G(4,1)=sin(conj(q(1)) + conj(q(2)) + conj(q(4)))*conj(g)*(conj(MT)*conj(lt) + 2*conj(Mf)*...
         conj(lt) + conj(Mfa)*conj(lt) + 2*conj(Mt)*conj(lt) + conj(Mua)*conj(lt) - conj(Mt)*conj(pMt));
  G(5,1)=-sin(conj(q(1)) + conj(q(3)) + conj(q(5)))*conj(Mt)*conj(g)*conj(pMt);
  G(6,1)=-conj(g)*(conj(Mfa)*(sin(conj(q(1)) + conj(q(6)))*conj(lua) + sin(conj(q(1)) +...
          conj(q(6)) - conj(q(7)))*conj(pMfa)) + sin(conj(q(1)) + conj(q(6)))*conj(Mua)*conj(pMua));
  G(7,1)=sin(conj(q(1)) + conj(q(6)) - conj(q(7)))*conj(Mfa)*conj(g)*conj(pMfa);

   D(1,1)=IT + 2*If + Ifa + 2*It + Iua + MT*lf^2 + MT*lt^2 + 2*Mf*lf^2 + Mfa*lf^2 + 2*Mf*lt^2 + 2*Mt*...
         lf^2 + Mfa*lt^2 + Mua*lf^2 + Mfa*lua^2 + 2*Mt*lt^2 + Mua*lt^2 + MT*pMT^2 + 2*Mf*pMf^2 + Mfa*pMfa^2 + 2*Mt*...
         pMt^2 + Mua*pMua^2 + Mfa*psh^2 + Mua*psh^2 - 2*Mt*lf^2*cos(q(2) - q(3)) - 2*Mf*lf*pMf - 2*Mt*lt*pMt - 2*MT*...
         lt*pMT*cos(q(2) + q(4)) - 2*Mfa*lt*psh*cos(q(2) + q(4)) - 2*Mua*lt*psh*cos(q(2) + q(4)) + 2*MT*lf*lt*...
         cos(q(4)) + 4*Mf*lf*lt*cos(q(4)) + 2*Mfa*lf*lt*cos(q(4)) + 2*Mt*lf*lt*cos(q(4)) + 2*Mua*lf*lt*cos(q(4)) - 2*MT*...
         lf*pMT*cos(q(2)) - 2*Mf*lt*pMf*cos(q(4)) + 2*Mt*lf*pMt*cos(q(5)) + 2*Mfa*lua*pMfa*cos(q(7)) - 2*Mfa*lf*...
         psh*cos(q(2)) - 2*Mua*lf*psh*cos(q(2)) + 2*Mfa*lua*psh*cos(q(6)) - 2*Mt*lf*pMt*cos(q(2) - q(3) - q(5)) +...
          2*Mua*pMua*psh*cos(q(6)) - 2*Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7)) - 2*Mfa*lf*lua*...
         cos(q(2) - q(6)) - 2*Mf*lf*pMf*cos(q(2) - q(3)) - 2*Mua*lf*pMua*cos(q(2) - q(6)) + 2*Mfa*pMfa*psh*cos(q(6) - q(7)) - 2*...
         Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5)) - 2*Mt*lf*lt*cos(q(2) - q(3) + q(4)) - 2*Mfa*lt*lua*cos(q(2) +...
          q(4) - q(6)) - 2*Mfa*lf*pMfa*cos(q(2) - q(6) + q(7)) - 2*Mf*lt*pMf*cos(q(2) - q(3) + q(4)) - 2*Mua*lt*pMua*cos(q(2) + q(4) - q(6));
  D(1,2)=If + It + MT*lf^2 + MT*lt^2 + 2*Mf*lf^2 + Mfa*lf^2 + 2*Mf*lt^2 + Mt*lf^2 + Mfa*lt^2 + Mua*...
         lf^2 + 2*Mt*lt^2 + Mua*lt^2 + Mf*pMf^2 + Mt*pMt^2 - Mt*lf^2*cos(q(2) - q(3)) - 2*Mf*lf*pMf - 2*Mt*lt*...
         pMt - MT*lt*pMT*cos(q(2) + q(4)) - Mfa*lt*psh*cos(q(2) + q(4)) - Mua*lt*psh*cos(q(2) + q(4)) + 2*MT*lf*lt*...
         cos(q(4)) + 4*Mf*lf*lt*cos(q(4)) + 2*Mfa*lf*lt*cos(q(4)) + 2*Mt*lf*lt*cos(q(4)) + 2*Mua*lf*lt*cos(q(4)) - MT*...
         lf*pMT*cos(q(2)) - 2*Mf*lt*pMf*cos(q(4)) - Mfa*lf*psh*cos(q(2)) - Mua*lf*psh*cos(q(2)) - Mt*lf*pMt*...
         cos(q(2) - q(3) - q(5)) - Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7)) - Mfa*lf*lua*cos(q(2) - q(6)) - Mf*lf*pMf*...
         cos(q(2) - q(3)) - Mua*lf*pMua*cos(q(2) - q(6)) - Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5)) - Mt*lf*lt*cos(q(2) - q(3) +...
          q(4)) - Mfa*lt*lua*cos(q(2) + q(4) - q(6)) - Mfa*lf*pMfa*cos(q(2) - q(6) + q(7)) - Mf*lt*pMf*cos(q(2) - q(3) +...
          q(4)) - Mua*lt*pMua*cos(q(2) + q(4) - q(6));
  D(1,3)=If + It + Mt*lf^2 + Mf*pMf^2 + Mt*pMt^2 - Mt*lf^2*cos(q(2) - q(3)) + 2*Mt*lf*pMt*...
         cos(q(5)) - Mt*lf*pMt*cos(q(2) - q(3) - q(5)) - Mf*lf*pMf*cos(q(2) - q(3)) - Mt*lt*pMt*cos(q(2) - q(3) +...
          q(4) - q(5)) - Mt*lf*lt*cos(q(2) - q(3) + q(4)) - Mf*lt*pMf*cos(q(2) - q(3) + q(4));
  D(1,4)=It + MT*lt^2 + 2*Mf*lt^2 + Mfa*lt^2 + 2*Mt*lt^2 + Mua*lt^2 + Mt*pMt^2 - 2*Mt*lt*pMt - MT*...
         lt*pMT*cos(q(2) + q(4)) - Mfa*lt*psh*cos(q(2) + q(4)) - Mua*lt*psh*cos(q(2) + q(4)) + MT*lf*lt*...
         cos(q(4)) + 2*Mf*lf*lt*cos(q(4)) + Mfa*lf*lt*cos(q(4)) + Mt*lf*lt*cos(q(4)) + Mua*lf*lt*cos(q(4)) - Mf*lt*pMf*...
         cos(q(4)) - Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7)) - Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5)) - Mt*lf*lt*...
         cos(q(2) - q(3) + q(4)) - Mfa*lt*lua*cos(q(2) + q(4) - q(6)) - Mf*lt*pMf*cos(q(2) - q(3) + q(4)) - Mua*lt*pMua*...
         cos(q(2) + q(4) - q(6));
  D(1,5)=It + Mt*pMt^2 + Mt*lf*pMt*cos(q(5)) - Mt*lf*pMt*cos(q(2) - q(3) - q(5)) - Mt*lt*pMt*...
         cos(q(2) - q(3) + q(4) - q(5));
  D(1,6)=Ifa + Iua + Mfa*lua^2 + Mfa*pMfa^2 + Mua*pMua^2 + 2*Mfa*lua*pMfa*cos(q(7)) + Mfa*lua*psh*...
         cos(q(6)) + Mua*pMua*psh*cos(q(6)) - Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7)) - Mfa*lf*lua*...
         cos(q(2) - q(6)) - Mua*lf*pMua*cos(q(2) - q(6)) + Mfa*pMfa*psh*cos(q(6) - q(7)) - Mfa*lt*lua*cos(q(2) + q(4) - q(6)) - Mfa*...
         lf*pMfa*cos(q(2) - q(6) + q(7)) - Mua*lt*pMua*cos(q(2) + q(4) - q(6));
  D(1,7)=Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7)) - Mfa*pMfa^2 - Mfa*lua*pMfa*...
         cos(q(7)) - Ifa - Mfa*pMfa*psh*cos(q(6) - q(7)) + Mfa*lf*pMfa*cos(q(2) - q(6) + q(7));
  D(2,1)=If + It + MT*lf^2 + MT*lt^2 + 2*Mf*lf^2 + Mfa*lf^2 + 2*Mf*lt^2 + Mt*lf^2 + Mfa*lt^2 + Mua*...
         lf^2 + 2*Mt*lt^2 + Mua*lt^2 + Mf*pMf^2 + Mt*pMt^2 - Mt*lf^2*cos(q(2) - q(3)) - 2*Mf*lf*pMf - 2*Mt*lt*...
         pMt - MT*lt*pMT*cos(q(2) + q(4)) - Mfa*lt*psh*cos(q(2) + q(4)) - Mua*lt*psh*cos(q(2) + q(4)) + 2*MT*lf*lt*...
         cos(q(4)) + 4*Mf*lf*lt*cos(q(4)) + 2*Mfa*lf*lt*cos(q(4)) + 2*Mt*lf*lt*cos(q(4)) + 2*Mua*lf*lt*cos(q(4)) - MT*...
         lf*pMT*cos(q(2)) - 2*Mf*lt*pMf*cos(q(4)) - Mfa*lf*psh*cos(q(2)) - Mua*lf*psh*cos(q(2)) - Mt*lf*pMt*...
         cos(q(2) - q(3) - q(5)) - Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7)) - Mfa*lf*lua*cos(q(2) - q(6)) - Mf*lf*pMf*...
         cos(q(2) - q(3)) - Mua*lf*pMua*cos(q(2) - q(6)) - Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5)) - Mt*lf*lt*cos(q(2) - q(3) +...
          q(4)) - Mfa*lt*lua*cos(q(2) + q(4) - q(6)) - Mfa*lf*pMfa*cos(q(2) - q(6) + q(7)) - Mf*lt*pMf*cos(q(2) - q(3) +...
          q(4)) - Mua*lt*pMua*cos(q(2) + q(4) - q(6));
  D(2,2)=If + It + MT*lf^2 + MT*lt^2 + 2*Mf*lf^2 + Mfa*lf^2 + 2*Mf*lt^2 + Mt*lf^2 + Mfa*lt^2 + Mua*...
         lf^2 + 2*Mt*lt^2 + Mua*lt^2 + Mf*pMf^2 + Mt*pMt^2 - 2*Mf*lf*pMf - 2*Mt*lt*pMt + 2*MT*lf*lt*cos(q(4)) + 4*...
         Mf*lf*lt*cos(q(4)) + 2*Mfa*lf*lt*cos(q(4)) + 2*Mt*lf*lt*cos(q(4)) + 2*Mua*lf*lt*cos(q(4)) - 2*Mf*lt*...
         pMf*cos(q(4));
  D(2,3)=- Mt*lf^2*cos(q(2) - q(3)) - Mt*lf*pMt*cos(q(2) - q(3) - q(5)) - Mf*lf*pMf*...
         cos(q(2) - q(3)) - Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5)) - Mt*lf*lt*cos(q(2) - q(3) + q(4)) - Mf*lt*pMf*...
         cos(q(2) - q(3) + q(4));
  D(2,4)=It + MT*lt^2 + 2*Mf*lt^2 + Mfa*lt^2 + 2*Mt*lt^2 + Mua*lt^2 + Mt*pMt^2 - 2*Mt*lt*pMt + MT*...
         lf*lt*cos(q(4)) + 2*Mf*lf*lt*cos(q(4)) + Mfa*lf*lt*cos(q(4)) + Mt*lf*lt*cos(q(4)) + Mua*lf*lt*...
         cos(q(4)) - Mf*lt*pMf*cos(q(4));
  D(2,5)=- Mt*lf*pMt*cos(q(2) - q(3) - q(5)) - Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5));
  D(2,6)=- Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7)) - Mfa*lf*lua*cos(q(2) - q(6)) - Mua*lf*pMua*...
         cos(q(2) - q(6)) - Mfa*lt*lua*cos(q(2) + q(4) - q(6)) - Mfa*lf*pMfa*cos(q(2) - q(6) + q(7)) - Mua*lt*pMua*cos(q(2) + q(4) - q(6));
  D(2,7)=Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7)) + Mfa*lf*pMfa*cos(q(2) - q(6) + q(7));
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
  D(3,6)=0;
  D(3,7)=0;
  D(4,1)=It + MT*lt^2 + 2*Mf*lt^2 + Mfa*lt^2 + 2*Mt*lt^2 + Mua*lt^2 + Mt*pMt^2 - 2*Mt*lt*pMt - MT*...
         lt*pMT*cos(q(2) + q(4)) - Mfa*lt*psh*cos(q(2) + q(4)) - Mua*lt*psh*cos(q(2) + q(4)) + MT*lf*lt*...
         cos(q(4)) + 2*Mf*lf*lt*cos(q(4)) + Mfa*lf*lt*cos(q(4)) + Mt*lf*lt*cos(q(4)) + Mua*lf*lt*cos(q(4)) - Mf*lt*pMf*...
         cos(q(4)) - Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7)) - Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5)) - Mt*lf*lt*...
         cos(q(2) - q(3) + q(4)) - Mfa*lt*lua*cos(q(2) + q(4) - q(6)) - Mf*lt*pMf*cos(q(2) - q(3) + q(4)) - Mua*lt*pMua*...
         cos(q(2) + q(4) - q(6));
  D(4,2)=It + MT*lt^2 + 2*Mf*lt^2 + Mfa*lt^2 + 2*Mt*lt^2 + Mua*lt^2 + Mt*pMt^2 - 2*Mt*lt*pMt + MT*...
         lf*lt*cos(q(4)) + 2*Mf*lf*lt*cos(q(4)) + Mfa*lf*lt*cos(q(4)) + Mt*lf*lt*cos(q(4)) + Mua*lf*lt*...
         cos(q(4)) - Mf*lt*pMf*cos(q(4));
  D(4,3)=- Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5)) - Mt*lf*lt*cos(q(2) - q(3) + q(4)) - Mf*lt*pMf*...
         cos(q(2) - q(3) + q(4));
  D(4,4)=It + MT*lt^2 + 2*Mf*lt^2 + Mfa*lt^2 + 2*Mt*lt^2 + Mua*lt^2 + Mt*pMt^2 - 2*Mt*lt*pMt;
  D(4,5)=-Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5));
  D(4,6)=- Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7)) - Mfa*lt*lua*cos(q(2) + q(4) - q(6)) - Mua*lt*...
         pMua*cos(q(2) + q(4) - q(6));
  D(4,7)=Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7));
  D(5,1)=It + Mt*pMt^2 + Mt*lf*pMt*cos(q(5)) - Mt*lf*pMt*cos(q(2) - q(3) - q(5)) - Mt*lt*pMt*...
         cos(q(2) - q(3) + q(4) - q(5));
  D(5,2)=- Mt*lf*pMt*cos(q(2) - q(3) - q(5)) - Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5));
  D(5,3)=It + Mt*pMt^2 + Mt*lf*pMt*cos(q(5));
  D(5,4)=-Mt*lt*pMt*cos(q(2) - q(3) + q(4) - q(5));
  D(5,5)=It + Mt*pMt^2;
  D(5,6)=0;
  D(5,7)=0;
  D(6,1)=Ifa + Iua + Mfa*lua^2 + Mfa*pMfa^2 + Mua*pMua^2 + 2*Mfa*lua*pMfa*cos(q(7)) + Mfa*lua*psh*...
         cos(q(6)) + Mua*pMua*psh*cos(q(6)) - Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7)) - Mfa*lf*lua*...
         cos(q(2) - q(6)) - Mua*lf*pMua*cos(q(2) - q(6)) + Mfa*pMfa*psh*cos(q(6) - q(7)) - Mfa*lt*lua*cos(q(2) + q(4) - q(6)) - Mfa*...
         lf*pMfa*cos(q(2) - q(6) + q(7)) - Mua*lt*pMua*cos(q(2) + q(4) - q(6));
  D(6,2)=- Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7)) - Mfa*lf*lua*cos(q(2) - q(6)) - Mua*lf*pMua*...
         cos(q(2) - q(6)) - Mfa*lt*lua*cos(q(2) + q(4) - q(6)) - Mfa*lf*pMfa*cos(q(2) - q(6) + q(7)) - Mua*lt*pMua*cos(q(2) + q(4) - q(6));
  D(6,3)=0;
  D(6,4)=- Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7)) - Mfa*lt*lua*cos(q(2) + q(4) - q(6)) - Mua*lt*...
         pMua*cos(q(2) + q(4) - q(6));
  D(6,5)=0;
  D(6,6)=Ifa + Iua + Mfa*lua^2 + Mfa*pMfa^2 + Mua*pMua^2 + 2*Mfa*lua*pMfa*cos(q(7));
  D(6,7)=- Ifa - Mfa*pMfa^2 - Mfa*lua*pMfa*cos(q(7));
  D(7,1)=Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7)) - Mfa*pMfa^2 - Mfa*lua*pMfa*...
         cos(q(7)) - Ifa - Mfa*pMfa*psh*cos(q(6) - q(7)) + Mfa*lf*pMfa*cos(q(2) - q(6) + q(7));
  D(7,2)=Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7)) + Mfa*lf*pMfa*cos(q(2) - q(6) + q(7));
  D(7,3)=0;
  D(7,4)=Mfa*lt*pMfa*cos(q(2) + q(4) - q(6) + q(7));
  D(7,5)=0;
  D(7,6)=- Ifa - Mfa*pMfa^2 - Mfa*lua*pMfa*cos(q(7));
  D(7,7)=Ifa + Mfa*pMfa^2;

   C(1,1)=Mt*dq(2)*lf^2*sin(q(2) - q(3)) - Mt*dq(3)*lf^2*sin(q(2) - q(3)) + Mf*dq(2)*lf*pMf*...
         sin(q(2) - q(3)) - Mf*dq(3)*lf*pMf*sin(q(2) - q(3)) + Mua*dq(2)*lf*pMua*sin(q(2) - q(6)) - Mua*dq(6)*lf*pMua*...
         sin(q(2) - q(6)) - Mfa*dq(6)*pMfa*psh*sin(q(6) - q(7)) + Mfa*dq(7)*pMfa*psh*sin(q(6) - q(7)) + Mt*dq(2)*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5)) - Mt*dq(3)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5)) + Mt*dq(4)*lt*pMt*sin(q(2) - q(3) +...
          q(4) - q(5)) - Mt*dq(5)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5)) + Mt*dq(2)*lf*lt*sin(q(2) - q(3) + q(4)) - Mt*dq(3)*lf*...
         lt*sin(q(2) - q(3) + q(4)) + Mt*dq(4)*lf*lt*sin(q(2) - q(3) + q(4)) + Mfa*dq(2)*lt*lua*sin(q(2) +...
          q(4) - q(6)) + Mfa*dq(4)*lt*lua*sin(q(2) + q(4) - q(6)) - Mfa*dq(6)*lt*lua*sin(q(2) + q(4) - q(6)) + Mfa*dq(2)*lf*...
         pMfa*sin(q(2) - q(6) + q(7)) - Mfa*dq(6)*lf*pMfa*sin(q(2) - q(6) + q(7)) + Mfa*dq(7)*lf*pMfa*...
         sin(q(2) - q(6) + q(7)) + Mf*dq(2)*lt*pMf*sin(q(2) - q(3) + q(4)) - Mf*dq(3)*lt*pMf*sin(q(2) - q(3) + q(4)) + Mf*...
         dq(4)*lt*pMf*sin(q(2) - q(3) + q(4)) + Mua*dq(2)*lt*pMua*sin(q(2) + q(4) - q(6)) + Mua*dq(4)*lt*pMua*...
         sin(q(2) + q(4) - q(6)) - Mua*dq(6)*lt*pMua*sin(q(2) + q(4) - q(6)) + MT*dq(2)*lt*pMT*sin(q(2) + q(4)) + MT*...
         dq(4)*lt*pMT*sin(q(2) + q(4)) + Mfa*dq(2)*lt*psh*sin(q(2) + q(4)) + Mfa*dq(4)*lt*psh*sin(q(2) + q(4)) +...
          Mua*dq(2)*lt*psh*sin(q(2) + q(4)) + Mua*dq(4)*lt*psh*sin(q(2) + q(4)) - MT*dq(4)*lf*lt*sin(q(4)) - 2*Mf*...
         dq(4)*lf*lt*sin(q(4)) - Mfa*dq(4)*lf*lt*sin(q(4)) - Mt*dq(4)*lf*lt*sin(q(4)) - Mua*dq(4)*lf*lt*sin(q(4)) +...
          MT*dq(2)*lf*pMT*sin(q(2)) + Mf*dq(4)*lt*pMf*sin(q(4)) - Mt*dq(5)*lf*pMt*sin(q(5)) - Mfa*dq(7)*lua*pMfa*...
         sin(q(7)) + Mfa*dq(2)*lf*psh*sin(q(2)) + Mua*dq(2)*lf*psh*sin(q(2)) - Mfa*dq(6)*lua*psh*sin(q(6)) + Mt*dq(2)*...
         lf*pMt*sin(q(2) - q(3) - q(5)) - Mt*dq(3)*lf*pMt*sin(q(2) - q(3) - q(5)) - Mt*dq(5)*lf*pMt*...
         sin(q(2) - q(3) - q(5)) - Mua*dq(6)*pMua*psh*sin(q(6)) + Mfa*dq(2)*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) + Mfa*dq(4)*lt*pMfa*...
         sin(q(2) + q(4) - q(6) + q(7)) - Mfa*dq(6)*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) + Mfa*dq(7)*lt*pMfa*...
         sin(q(2) + q(4) - q(6) + q(7)) + Mfa*dq(2)*lf*lua*sin(q(2) - q(6)) - Mfa*dq(6)*lf*lua*sin(q(2) - q(6));
  C(1,2)=Mt*dq(1)*lf^2*sin(q(2) - q(3)) + Mt*dq(2)*lf^2*sin(q(2) - q(3)) + Mf*dq(1)*lf*pMf*...
         sin(q(2) - q(3)) + Mf*dq(2)*lf*pMf*sin(q(2) - q(3)) + Mua*dq(1)*lf*pMua*sin(q(2) - q(6)) + Mua*dq(2)*lf*pMua*...
         sin(q(2) - q(6)) + Mt*dq(1)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5)) + Mt*dq(2)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5)) +...
          Mt*dq(4)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5)) + Mt*dq(1)*lf*lt*sin(q(2) - q(3) + q(4)) + Mt*dq(2)*lf*...
         lt*sin(q(2) - q(3) + q(4)) + Mt*dq(4)*lf*lt*sin(q(2) - q(3) + q(4)) + Mfa*dq(1)*lt*lua*sin(q(2) +...
          q(4) - q(6)) + Mfa*dq(2)*lt*lua*sin(q(2) + q(4) - q(6)) + Mfa*dq(4)*lt*lua*sin(q(2) + q(4) - q(6)) + Mfa*dq(1)*lf*...
         pMfa*sin(q(2) - q(6) + q(7)) + Mfa*dq(2)*lf*pMfa*sin(q(2) - q(6) + q(7)) + Mf*dq(1)*lt*pMf*...
         sin(q(2) - q(3) + q(4)) + Mf*dq(2)*lt*pMf*sin(q(2) - q(3) + q(4)) + Mf*dq(4)*lt*pMf*sin(q(2) - q(3) + q(4)) + Mua*...
         dq(1)*lt*pMua*sin(q(2) + q(4) - q(6)) + Mua*dq(2)*lt*pMua*sin(q(2) + q(4) - q(6)) + Mua*dq(4)*lt*pMua*...
         sin(q(2) + q(4) - q(6)) + MT*dq(1)*lt*pMT*sin(q(2) + q(4)) + MT*dq(2)*lt*pMT*sin(q(2) + q(4)) + MT*dq(4)*lt*...
         pMT*sin(q(2) + q(4)) + Mfa*dq(1)*lt*psh*sin(q(2) + q(4)) + Mfa*dq(2)*lt*psh*sin(q(2) + q(4)) + Mfa*...
         dq(4)*lt*psh*sin(q(2) + q(4)) + Mua*dq(1)*lt*psh*sin(q(2) + q(4)) + Mua*dq(2)*lt*psh*sin(q(2) + q(4)) +...
          Mua*dq(4)*lt*psh*sin(q(2) + q(4)) - MT*dq(4)*lf*lt*sin(q(4)) - 2*Mf*dq(4)*lf*lt*sin(q(4)) - Mfa*dq(4)*...
         lf*lt*sin(q(4)) - Mt*dq(4)*lf*lt*sin(q(4)) - Mua*dq(4)*lf*lt*sin(q(4)) + MT*dq(1)*lf*pMT*sin(q(2)) +...
          MT*dq(2)*lf*pMT*sin(q(2)) + Mf*dq(4)*lt*pMf*sin(q(4)) + Mfa*dq(1)*lf*psh*sin(q(2)) + Mfa*dq(2)*lf*psh*...
         sin(q(2)) + Mua*dq(1)*lf*psh*sin(q(2)) + Mua*dq(2)*lf*psh*sin(q(2)) + Mt*dq(1)*lf*pMt*sin(q(2) - q(3) - q(5)) +...
          Mt*dq(2)*lf*pMt*sin(q(2) - q(3) - q(5)) + Mfa*dq(1)*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) + Mfa*dq(2)*...
         lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) + Mfa*dq(4)*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) + Mfa*dq(1)*...
         lf*lua*sin(q(2) - q(6)) + Mfa*dq(2)*lf*lua*sin(q(2) - q(6));
  C(1,3)=- dq(1)*(Mt*lf^2*sin(q(2) - q(3)) + Mt*lf*lt*sin(q(2) - q(3) + q(4)) + Mf*lt*pMf*...
         sin(q(2) - q(3) + q(4)) + Mt*lf*pMt*sin(q(2) - q(3) - q(5)) + Mf*lf*pMf*sin(q(2) - q(3)) + Mt*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5))) - dq(3)*(Mt*lf^2*sin(q(2) - q(3)) + Mt*lf*lt*sin(q(2) - q(3) + q(4)) + Mf*lt*pMf*...
         sin(q(2) - q(3) + q(4)) + Mt*lf*pMt*sin(q(2) - q(3) - q(5)) + Mf*lf*pMf*sin(q(2) - q(3)) + Mt*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5))) - Mt*dq(5)*pMt*(lt*sin(q(2) - q(3) + q(4) - q(5)) + lf*sin(q(5)) + lf*sin(q(2) - q(3) - q(5)));
  C(1,4)=lt*(dq(1) + dq(2) + dq(4))*(Mf*pMf*sin(q(4)) + Mfa*pMfa*sin(q(2) + q(4) - q(6) + q(7)) +...
          Mt*pMt*sin(q(2) - q(3) + q(4) - q(5)) + Mt*lf*sin(q(2) - q(3) + q(4)) + Mfa*lua*sin(q(2) +...
          q(4) - q(6)) + Mf*pMf*sin(q(2) - q(3) + q(4)) + Mua*pMua*sin(q(2) + q(4) - q(6)) + MT*pMT*sin(q(2) + q(4)) + Mfa*...
         psh*sin(q(2) + q(4)) + Mua*psh*sin(q(2) + q(4)) - MT*lf*sin(q(4)) - 2*Mf*lf*sin(q(4)) - Mfa*lf*...
         sin(q(4)) - Mt*lf*sin(q(4)) - Mua*lf*sin(q(4)));
  C(1,5)=-Mt*pMt*(dq(1) + dq(3) + dq(5))*(lt*sin(q(2) - q(3) + q(4) - q(5)) + lf*sin(q(5)) + lf*sin(q(2) - q(3) - q(5)));
  C(1,6)=Mfa*dq(7)*pMfa*(lf*sin(q(2) - q(6) + q(7)) - lua*sin(q(7)) + lt*sin(q(2) + q(4) - q(6) +...
          q(7)) + psh*sin(q(6) - q(7))) - dq(6)*(Mfa*lt*lua*sin(q(2) + q(4) - q(6)) + Mfa*lf*pMfa*sin(q(2) - q(6) +...
          q(7)) + Mua*lt*pMua*sin(q(2) + q(4) - q(6)) + Mfa*lua*psh*sin(q(6)) + Mua*pMua*psh*sin(q(6)) + Mfa*lt*pMfa*...
         sin(q(2) + q(4) - q(6) + q(7)) + Mfa*lf*lua*sin(q(2) - q(6)) + Mua*lf*pMua*sin(q(2) - q(6)) + Mfa*pMfa*psh*...
         sin(q(6) - q(7))) - dq(1)*(Mfa*lt*lua*sin(q(2) + q(4) - q(6)) + Mfa*lf*pMfa*sin(q(2) - q(6) + q(7)) + Mua*lt*pMua*sin(q(2) +...
          q(4) - q(6)) + Mfa*lua*psh*sin(q(6)) + Mua*pMua*psh*sin(q(6)) + Mfa*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) + Mfa*...
         lf*lua*sin(q(2) - q(6)) + Mua*lf*pMua*sin(q(2) - q(6)) + Mfa*pMfa*psh*sin(q(6) - q(7)));
  C(1,7)=Mfa*pMfa*(dq(1) + dq(6) - dq(7))*(lf*sin(q(2) - q(6) + q(7)) - lua*sin(q(7)) + lt*...
         sin(q(2) + q(4) - q(6) + q(7)) + psh*sin(q(6) - q(7)));
  C(2,1)=Mfa*dq(7)*lf*pMfa*sin(q(2) - q(6) + q(7)) - Mt*dq(3)*lf^2*sin(q(2) - q(3)) - Mf*dq(1)*lf*...
         pMf*sin(q(2) - q(3)) - Mf*dq(3)*lf*pMf*sin(q(2) - q(3)) - Mua*dq(1)*lf*pMua*sin(q(2) - q(6)) - Mua*...
         dq(6)*lf*pMua*sin(q(2) - q(6)) - Mt*dq(1)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5)) - Mt*dq(3)*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5)) - Mt*dq(5)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5)) - Mt*dq(1)*lf*lt*sin(q(2) - q(3) +...
          q(4)) - Mt*dq(3)*lf*lt*sin(q(2) - q(3) + q(4)) - Mfa*dq(1)*lt*lua*sin(q(2) + q(4) - q(6)) - Mfa*dq(6)*lt*lua*...
         sin(q(2) + q(4) - q(6)) - Mfa*dq(1)*lf*pMfa*sin(q(2) - q(6) + q(7)) - Mfa*dq(6)*lf*pMfa*sin(q(2) - q(6) +...
          q(7)) - Mt*dq(1)*lf^2*sin(q(2) - q(3)) - Mf*dq(1)*lt*pMf*sin(q(2) - q(3) + q(4)) - Mf*dq(3)*lt*pMf*...
         sin(q(2) - q(3) + q(4)) - Mua*dq(1)*lt*pMua*sin(q(2) + q(4) - q(6)) - Mua*dq(6)*lt*pMua*sin(q(2) + q(4) - q(6)) - MT*...
         dq(1)*lt*pMT*sin(q(2) + q(4)) - Mfa*dq(1)*lt*psh*sin(q(2) + q(4)) - Mua*dq(1)*lt*psh*sin(q(2) +...
          q(4)) - MT*dq(4)*lf*lt*sin(q(4)) - 2*Mf*dq(4)*lf*lt*sin(q(4)) - Mfa*dq(4)*lf*lt*sin(q(4)) - Mt*dq(4)*lf*lt*...
         sin(q(4)) - Mua*dq(4)*lf*lt*sin(q(4)) - MT*dq(1)*lf*pMT*sin(q(2)) + Mf*dq(4)*lt*pMf*sin(q(4)) - Mfa*dq(1)*lf*psh*...
         sin(q(2)) - Mua*dq(1)*lf*psh*sin(q(2)) - Mt*dq(1)*lf*pMt*sin(q(2) - q(3) - q(5)) - Mt*dq(3)*lf*pMt*...
         sin(q(2) - q(3) - q(5)) - Mt*dq(5)*lf*pMt*sin(q(2) - q(3) - q(5)) - Mfa*dq(1)*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) - Mfa*dq(6)*...
         lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) + Mfa*dq(7)*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) - Mfa*dq(1)*...
         lf*lua*sin(q(2) - q(6)) - Mfa*dq(6)*lf*lua*sin(q(2) - q(6));
  C(2,2)=-dq(4)*lt*sin(q(4))*(MT*lf + 2*Mf*lf + Mfa*lf + Mt*lf + Mua*lf - Mf*pMf);
  C(2,3)=- dq(1)*(Mt*lf^2*sin(q(2) - q(3)) + Mt*lf*lt*sin(q(2) - q(3) + q(4)) + Mf*lt*pMf*...
         sin(q(2) - q(3) + q(4)) + Mt*lf*pMt*sin(q(2) - q(3) - q(5)) + Mf*lf*pMf*sin(q(2) - q(3)) + Mt*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5))) - dq(3)*(Mt*lf^2*sin(q(2) - q(3)) + Mt*lf*lt*sin(q(2) - q(3) + q(4)) + Mf*lt*pMf*...
         sin(q(2) - q(3) + q(4)) + Mt*lf*pMt*sin(q(2) - q(3) - q(5)) + Mf*lf*pMf*sin(q(2) - q(3)) + Mt*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5))) - Mt*dq(5)*pMt*(lt*sin(q(2) - q(3) + q(4) - q(5)) + lf*sin(q(2) - q(3) - q(5)));
  C(2,4)=-lt*sin(q(4))*(dq(1) + dq(2) + dq(4))*(MT*lf + 2*Mf*lf + Mfa*lf + Mt*lf + Mua*lf - Mf*pMf);
  C(2,5)=-Mt*pMt*(lt*sin(q(2) - q(3) + q(4) - q(5)) + lf*sin(q(2) - q(3) - q(5)))*(dq(1) + dq(3) + dq(5));
  C(2,6)=Mfa*dq(7)*pMfa*(lf*sin(q(2) - q(6) + q(7)) + lt*sin(q(2) + q(4) - q(6) + q(7))) - dq(6)*...
         (Mfa*lt*lua*sin(q(2) + q(4) - q(6)) + Mfa*lf*pMfa*sin(q(2) - q(6) + q(7)) + Mua*lt*pMua*sin(q(2) +...
          q(4) - q(6)) + Mfa*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) + Mfa*lf*lua*sin(q(2) - q(6)) + Mua*lf*pMua*...
         sin(q(2) - q(6))) - dq(1)*(Mfa*lt*lua*sin(q(2) + q(4) - q(6)) + Mfa*lf*pMfa*sin(q(2) - q(6) + q(7)) + Mua*lt*pMua*sin(q(2) +...
          q(4) - q(6)) + Mfa*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) + Mfa*lf*lua*sin(q(2) - q(6)) + Mua*lf*pMua*sin(q(2) - q(6)));
  C(2,7)=Mfa*pMfa*(lf*sin(q(2) - q(6) + q(7)) + lt*sin(q(2) + q(4) - q(6) + q(7)))*(dq(1) + dq(6) - dq(7));
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
  C(3,6)=0;
  C(3,7)=0;
  C(4,1)=MT*dq(1)*lf*lt*sin(q(4)) - Mt*dq(3)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5)) - Mt*dq(5)*lt*...
         pMt*sin(q(2) - q(3) + q(4) - q(5)) - Mt*dq(1)*lf*lt*sin(q(2) - q(3) + q(4)) - Mt*dq(3)*lf*lt*...
         sin(q(2) - q(3) + q(4)) - Mfa*dq(1)*lt*lua*sin(q(2) + q(4) - q(6)) - Mfa*dq(6)*lt*lua*sin(q(2) + q(4) - q(6)) - Mf*...
         dq(1)*lt*pMf*sin(q(2) - q(3) + q(4)) - Mf*dq(3)*lt*pMf*sin(q(2) - q(3) + q(4)) - Mua*dq(1)*lt*pMua*...
         sin(q(2) + q(4) - q(6)) - Mua*dq(6)*lt*pMua*sin(q(2) + q(4) - q(6)) - MT*dq(1)*lt*pMT*sin(q(2) + q(4)) - Mfa*...
         dq(1)*lt*psh*sin(q(2) + q(4)) - Mua*dq(1)*lt*psh*sin(q(2) + q(4)) - Mt*dq(1)*lt*pMt*sin(q(2) - q(3) +...
          q(4) - q(5)) + MT*dq(2)*lf*lt*sin(q(4)) + 2*Mf*dq(1)*lf*lt*sin(q(4)) + 2*Mf*dq(2)*lf*lt*sin(q(4)) + Mfa*dq(1)*lf*...
         lt*sin(q(4)) + Mfa*dq(2)*lf*lt*sin(q(4)) + Mt*dq(1)*lf*lt*sin(q(4)) + Mt*dq(2)*lf*lt*sin(q(4)) + Mua*...
         dq(1)*lf*lt*sin(q(4)) + Mua*dq(2)*lf*lt*sin(q(4)) - Mf*dq(1)*lt*pMf*sin(q(4)) - Mf*dq(2)*lt*pMf*...
         sin(q(4)) - Mfa*dq(1)*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) - Mfa*dq(6)*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) +...
          Mfa*dq(7)*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7));
  C(4,2)=lt*sin(q(4))*(dq(1) + dq(2))*(MT*lf + 2*Mf*lf + Mfa*lf + Mt*lf + Mua*lf - Mf*pMf);
  C(4,3)=- dq(1)*lt*(Mt*pMt*sin(q(2) - q(3) + q(4) - q(5)) + Mt*lf*sin(q(2) - q(3) + q(4)) + Mf*pMf*...
         sin(q(2) - q(3) + q(4))) - dq(3)*lt*(Mt*pMt*sin(q(2) - q(3) + q(4) - q(5)) + Mt*lf*sin(q(2) - q(3) + q(4)) + Mf*pMf*...
         sin(q(2) - q(3) + q(4))) - Mt*dq(5)*lt*pMt*sin(q(2) - q(3) + q(4) - q(5));
  C(4,4)=0;
  C(4,5)=-Mt*lt*pMt*sin(q(2) - q(3) + q(4) - q(5))*(dq(1) + dq(3) + dq(5));
  C(4,6)=Mfa*dq(7)*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) - dq(6)*lt*(Mfa*pMfa*sin(q(2) +...
          q(4) - q(6) + q(7)) + Mfa*lua*sin(q(2) + q(4) - q(6)) + Mua*pMua*sin(q(2) + q(4) - q(6))) - dq(1)*lt*(Mfa*pMfa*...
         sin(q(2) + q(4) - q(6) + q(7)) + Mfa*lua*sin(q(2) + q(4) - q(6)) + Mua*pMua*sin(q(2) + q(4) - q(6)));
  C(4,7)=Mfa*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7))*(dq(1) + dq(6) - dq(7));
  C(5,1)=Mt*pMt*(dq(1)*lf*sin(q(5)) + dq(3)*lf*sin(q(5)) + dq(1)*lf*sin(q(2) - q(3) - q(5)) + dq(2)*...
         lf*sin(q(2) - q(3) - q(5)) + dq(1)*lt*sin(q(2) - q(3) + q(4) - q(5)) + dq(2)*lt*sin(q(2) - q(3) +...
          q(4) - q(5)) + dq(4)*lt*sin(q(2) - q(3) + q(4) - q(5)));
  C(5,2)=dq(1)*(Mt*lf*pMt*sin(q(2) - q(3) - q(5)) + Mt*lt*pMt*sin(q(2) - q(3) + q(4) - q(5))) +...
          dq(2)*(Mt*lf*pMt*sin(q(2) - q(3) - q(5)) + Mt*lt*pMt*sin(q(2) - q(3) + q(4) - q(5))) + Mt*dq(4)*lt*pMt*...
         sin(q(2) - q(3) + q(4) - q(5));
  C(5,3)=Mt*lf*pMt*sin(q(5))*(dq(1) + dq(3));
  C(5,4)=Mt*lt*pMt*sin(q(2) - q(3) + q(4) - q(5))*(dq(1) + dq(2) + dq(4));
  C(5,5)=0;
  C(5,6)=0;
  C(5,7)=0;
  C(6,1)=Mua*dq(1)*lf*pMua*sin(q(2) - q(6)) + Mua*dq(2)*lf*pMua*sin(q(2) - q(6)) + Mfa*dq(1)*pMfa*...
         psh*sin(q(6) - q(7)) + Mfa*dq(1)*lt*lua*sin(q(2) + q(4) - q(6)) + Mfa*dq(2)*lt*lua*sin(q(2) +...
          q(4) - q(6)) + Mfa*dq(4)*lt*lua*sin(q(2) + q(4) - q(6)) + Mfa*dq(1)*lf*pMfa*sin(q(2) - q(6) + q(7)) + Mfa*dq(2)*...
         lf*pMfa*sin(q(2) - q(6) + q(7)) + Mua*dq(1)*lt*pMua*sin(q(2) + q(4) - q(6)) + Mua*dq(2)*lt*pMua*...
         sin(q(2) + q(4) - q(6)) + Mua*dq(4)*lt*pMua*sin(q(2) + q(4) - q(6)) - Mfa*dq(7)*lua*pMfa*sin(q(7)) + Mfa*...
         dq(1)*lua*psh*sin(q(6)) + Mua*dq(1)*pMua*psh*sin(q(6)) + Mfa*dq(1)*lt*pMfa*sin(q(2) + q(4) - q(6) +...
          q(7)) + Mfa*dq(2)*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) + Mfa*dq(4)*lt*pMfa*sin(q(2) + q(4) - q(6) +...
          q(7)) + Mfa*dq(1)*lf*lua*sin(q(2) - q(6)) + Mfa*dq(2)*lf*lua*sin(q(2) - q(6));
  C(6,2)=dq(1)*(Mfa*lt*lua*sin(q(2) + q(4) - q(6)) + Mfa*lf*pMfa*sin(q(2) - q(6) + q(7)) + Mua*lt*...
         pMua*sin(q(2) + q(4) - q(6)) + Mfa*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) + Mfa*lf*lua*sin(q(2) - q(6)) +...
          Mua*lf*pMua*sin(q(2) - q(6))) + dq(2)*(Mfa*lt*lua*sin(q(2) + q(4) - q(6)) + Mfa*lf*pMfa*...
         sin(q(2) - q(6) + q(7)) + Mua*lt*pMua*sin(q(2) + q(4) - q(6)) + Mfa*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7)) + Mfa*lf*...
         lua*sin(q(2) - q(6)) + Mua*lf*pMua*sin(q(2) - q(6))) + dq(4)*lt*(Mfa*pMfa*sin(q(2) + q(4) - q(6) +...
          q(7)) + Mfa*lua*sin(q(2) + q(4) - q(6)) + Mua*pMua*sin(q(2) + q(4) - q(6)));
  C(6,3)=0;
  C(6,4)=lt*(dq(1) + dq(2) + dq(4))*(Mfa*pMfa*sin(q(2) + q(4) - q(6) + q(7)) + Mfa*lua*sin(q(2) +...
          q(4) - q(6)) + Mua*pMua*sin(q(2) + q(4) - q(6)));
  C(6,5)=0;
  C(6,6)=-Mfa*dq(7)*lua*pMfa*sin(q(7));
  C(6,7)=-Mfa*lua*pMfa*sin(q(7))*(dq(1) + dq(6) - dq(7));
  C(7,1)=-Mfa*pMfa*(dq(1)*lt*sin(q(2) + q(4) - q(6) + q(7)) - dq(6)*lua*sin(q(7)) - dq(1)*lua*...
         sin(q(7)) + dq(2)*lt*sin(q(2) + q(4) - q(6) + q(7)) + dq(4)*lt*sin(q(2) + q(4) - q(6) + q(7)) + dq(1)*psh*...
         sin(q(6) - q(7)) + dq(1)*lf*sin(q(2) - q(6) + q(7)) + dq(2)*lf*sin(q(2) - q(6) + q(7)));
  C(7,2)=- dq(1)*(Mfa*lf*pMfa*sin(q(2) - q(6) + q(7)) + Mfa*lt*pMfa*sin(q(2) + q(4) - q(6) +...
          q(7))) - dq(2)*(Mfa*lf*pMfa*sin(q(2) - q(6) + q(7)) + Mfa*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7))) - Mfa*dq(4)*lt*...
         pMfa*sin(q(2) + q(4) - q(6) + q(7));
  C(7,3)=0;
  C(7,4)=-Mfa*lt*pMfa*sin(q(2) + q(4) - q(6) + q(7))*(dq(1) + dq(2) + dq(4));
  C(7,5)=0;
  C(7,6)=Mfa*lua*pMfa*sin(q(7))*(dq(1) + dq(6));
  C(7,7)=0;

 