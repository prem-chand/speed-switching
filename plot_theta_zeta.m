clear all;

fps = load('./Fixed_Point_Lib/fp_F=0/fps_lib_fine.mat');
data = fps.data;

fix=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0.mat');
ival=fix.x;

jacob=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0_time_beta_jacobian.mat');
G = jacob.G;

for i = 1:1:size(data,1)
    c = i/size(data,1);
    [theta,zeta] = compute_theta_zeta(data(i,2),ival,G,data(i,3:12)');
    speed = data(i,1)*ones([1,length(theta)]);
    h1 = plot3(theta,zeta,speed,'color',[c 1-c 0]);
    set(h1,'linewidth',3)
    hold on;
    h2 = plot3([theta(1) theta(end)],[zeta(1) zeta(end)],[speed(1) speed(end)],'color',[c 1-c 0]);
    set(h2,'linewidth',3)
%     if i==size(data,1)
%         plot(theta,zeta,'LineWidth',3,'Color',[1 0 0]);
%         hold on;
%         plot([theta(1) theta(end)],[zeta(1) zeta(end)],'LineWidth',3,'Color',[1 0 0]);
%         %     elseif i==1
%         %         plot(theta,zeta,'LineWidth',3,'Color',[0 0 1]);
%         %         hold on;
%         %         plot([theta(1) theta(end)],[zeta(1) zeta(end)],'--','LineWidth',3,'Color',[0 0 1]);
%     else
%         plot(theta,zeta,'LineWidth',2,'Color',[0.8 0.8 0.8]);
%         hold on;
%         plot([theta(1) theta(end)],[zeta(1) zeta(end)],'LineWidth',2,'Color',[0.8 0.8 0.8]);
%     end
    i
    clear theta zeta speed
end
% [theta,zeta] = compute_theta_zeta(data(1,2),ival,G,data(1,3:12)');
% plot(theta,zeta,'LineWidth',3,'Color',[0 0 1]);
% hold on;
% plot([theta(1) theta(end)],[zeta(1) zeta(end)],'LineWidth',3,'Color',[0 0 1]);
hold off;