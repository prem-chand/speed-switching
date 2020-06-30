function dwell_time = compute_dwell_time(epsilon,fp1,fp2,delta_z)
%COMPUTE_FP_SPEED_BETA Summary of this function goes here
%   Detailed explanation goes here
dwell_time = ceil(0.5*(log((abs(fp1-fp2)/epsilon)+1))/(log(1/delta_z)));
end