function Ff=fcn_out_force(q,u)
        dq=[q(6); q(7); q(8); q(9); q(10); 0;0 ];

        [D_e,C_e,G_e,JF,dJF] = fcn_extended_DCG(q);
        B_e=[zeros(1,4);eye(4);zeros(2,4)];

        invD_e = inv(D_e);
        Ff = ( inv(JF*invD_e*JF') * ...
            ( (JF*invD_e*C_e - dJF)*dq + JF*invD_e*G_e - ...
            JF*invD_e*B_e*u ))';
        Ff = Ff';
end
        