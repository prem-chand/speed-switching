function Phead=fcn_position_head(q)

lt=0.4;lf=0.4;lT=0.63;

Phead=[-lt*sin(pi-q(1)-q(2)-q(4))-lf*sin(pi-q(1)-q(2))+lT*sin(q(1)) lt*cos(pi-q(1)-q(2)-q(4))+lf*cos(pi-q(1)-q(2))+lT*cos(q(1))];
end