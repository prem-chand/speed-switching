clear all

fps = load('./Fixed_Point_Lib/fp_F=0/fps_lib_coarse.mat');
data = fps.data;

fix=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0.mat');
ival=fix.x;

jacob=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0_time_beta_jacobian.mat');
G = jacob.G;

edge = load('./Fixed_Point_Lib/fp_F=0/edge_coarse_PD.mat');
edge_matrix_coarse = edge.edge_matrix;

% index = [1 2 50 65 79];
index = 1:1:size(data,1);

[X,Y] = meshgrid(index,index);

X = reshape(X,[length(index)*length(index) 1]);
Y = reshape(Y,[length(index)*length(index) 1]);
% edge_matrix = eye(length(index));

speed_ctrl = data(:,2);
fp_start = data(:,3:12);
orbit_speed = data(:,1);

tic
parfor_progress(length(X));
parfor i=1:1:length(X)
    if X(i)~=Y(i) && edge_matrix_coarse(X(i),Y(i))==0
        check = compute_fp_trans_u_cone(speed_ctrl(Y(i)),ival,G,fp_start(X(i),:)');
%         speed = compute_speed_check_convergence(speed_ctrl(Y(i)),ival,G,fp_start(X(i),:)');
%         if (abs(orbit_speed(Y(i))-speed)<=0.025)
%             check = 1;
%         else check = 0;
%         end
        edge_matrix(i) = check;
    else
        edge_matrix(i)=1;
    end
    parfor_progress;
end
toc
parfor_progress(0);
edge_matrix = reshape(edge_matrix,[length(index) length(index)]);
edge_matrix = edge_matrix';
toc

save('./Fixed_Point_Lib/fp_F=0/edge_coarse_CLF_108.mat','edge_matrix');