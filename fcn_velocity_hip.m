function Ph=fcn_velocity_hip(q)
lt=0.4;lf=0.4;
Ph=-lt*(-q(6)-q(7)-q(9))*cos(pi-q(1)-q(2)-q(4))-(-q(6)-q(7))*lf*cos(pi-q(1)-q(2));
end