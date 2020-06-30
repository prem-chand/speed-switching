function [x,Fground]=impact_map(q)

[D_e,E_e]=fcn_extended_D_E(q);

dqe_minus=[q(6) q(7) q(8) q(9) q(10) 0 0]';
dqe_plus=[D_e -E_e';E_e zeros(2,2)]\[D_e*dqe_minus;0;0];
trans=inv([D_e -E_e';E_e zeros(2,2)]);
x=[q(1);q(3);q(2);q(5);q(4);dqe_plus(1);dqe_plus(3);dqe_plus(2);dqe_plus(5);dqe_plus(4)];
Fground=[dqe_plus(8) dqe_plus(9)];

end