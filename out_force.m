function F=out_force(x,u)

        [D_e,C_e,G_e,JF,dJF] = fcn_extended_DCG(x);

        invD_e = inv(D_e);

        F = ( inv(JF*invD_e*JF') * ...
            ( (JF*invD_e*C_e - dJF)*dq_e + JF*invD_e*G_e - ...
            JF*invD_e*B_e*u' ))';
        