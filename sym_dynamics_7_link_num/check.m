[G,D,C]=fcn_dynamics_matrices2(rand([1 14]));
Dl=D(1:5,1:5);
Dm=D(6:7,6:7);
Dml=D(1:5,6:7);
check_det=(eye(5,5)-inv(Dl)*Dml*inv(Dm)*transpose(Dml))
det(check_det)