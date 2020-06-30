function Pk=fcn_position_stanceknee(q)
lt=0.4;lf=0.4;
Pk=[-lt*sin(pi-q(1)-q(2)-q(4)) lt*cos(pi-q(1)-q(2)-q(4))];
end
