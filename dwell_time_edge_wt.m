fps = load('./Fixed_Point_Lib/fp_F=0/fps_lib_coarse.mat');
data = fps.data;

edge = load('./Fixed_Point_Lib/fp_F=0/edge_coarse_CLF_108.mat');
edge_matrix = edge.edge_matrix;

epsilon = 2;
delta_z = 0.9159;

for i = 1:1:size(edge_matrix,1)
    for j = 1:1:size(edge_matrix,2)
        if (edge_matrix(i,j) == 1)
            edge_matrix(i,j) = compute_dwell_time(epsilon,data(i,13),data(j,13),delta_z);
        end
    end
end
save('./Fixed_Point_Lib/fp_F=0/edge_coarse_CLF_108_wts.mat','edge_matrix');