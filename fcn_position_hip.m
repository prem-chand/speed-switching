function Ph=fcn_position_hip(q)
lt=0.4;lf=0.4;
Ph=[-lt*sin(pi-q(1)-q(2)-q(4))-lf*sin(pi-q(1)-q(2)) lt*cos(pi-q(1)-q(2)-q(4))+lf*cos(pi-q(1)-q(2))];
end

