clear all;
close all;
fix=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0.mat');
ival=fix.x;
XX=ival(1:10)';
jacob=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0_time_beta_jacobian.mat');
G = jacob.G;
speed_ctrl = [(-31.5:0.5:7.0)*1e-3 7.4e-3];
parfor_progress(length(speed_ctrl));
parfor i = 1:1:length(speed_ctrl)
    [avg_speed(i),fp(:,i)] = compute_fp_speed_beta(speed_ctrl(i),ival,G);
    [zeta_minus_fp(i),zeta_minus_lb(i)] = compute_zeta_limit(speed_ctrl(i),ival,G,fp(:,i));
    parfor_progress;
end
parfor_progress(0);
data = [avg_speed' speed_ctrl' fp' zeta_minus_fp' zeta_minus_lb'];
save('./Fixed_Point_Lib/fp_F=0/fps_lib_coarse_CLF.mat','data');