clear all;
close all;
clc;

%% Initialization
M=6;
N = 8;
te=0;
tplot=0;timpact=0;
ntouch=0;ptouch=0;
accel(1)=0;
simout.t=[];simout.stk=[];simout.h=[]; simout.swk=[]; simout.swf=[]; simout.head=[];
tF=[0 500];
Fext=[]; % No external force here
start = 79;
path = start;
warning('off','all');

global theta_minus theta_plus betta M_bez

fix=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0.mat');
ival=fix.x; % Encodes controller parameters for the "base" limit cycle

edge = load('./Fixed_Point_Lib/fp_F=0/edge_coarse_CLF_108_wts.mat');
edge_matrix_coarse = edge.edge_matrix;
gr = digraph(edge_matrix_coarse); % Generates a digraph for allowable gait switches
bins = conncomp(gr);

fps = load('./Fixed_Point_Lib/fp_F=0/fps_lib_coarse.mat');
data = fps.data; % Library of fixed points; see readme for data structure
XX = data(start,3:12)'; % Start is the index of the starting fp in data
desired_speed = [0.42 0.81]; % Give sequence of speeds desired to be achieved in a row here
for i = 1:1:length(desired_speed)
    [~,desired_point(i)] = min(abs(data(:,1)-desired_speed(i))); %  Find the node on the graph closest to th desired speed
    [path_intermediate,cost(i)] = shortestpath(gr,path(end),desired_point(i)); % Find the shortest path based on edge weight
    path = [path path_intermediate(2:end)];
end
sequence = data(path,1)' % Output the sequence of speed switches required to achieve the desired speeds
k = start;
zeta = data(start,13); % Output used to compute the dwell-time (can be interpreted as a zero dynamic state)

jacob=open('fixedpointforfivelink/Switching_Control/fixed_point_F=0_time_beta_jacobian.mat');
G = jacob.G; % Average speed jacobian w.r.t control parameters beta, see Section V of Veer et al CDC 2017

M_bez=6;
c=[1 1 0 0.5 0];
dwell_time = 0;
last_step = 0;
count = 2;
delta_z = 0.9159;
epsilon = 2;
i = 1;
count2 = 1;
t_history = [];
u_history = [];

%%  Execution
while count~=-1
    clear s theta_minus theta_plus betta event
    if i >= last_step + dwell_time
        if count>length(path)
            break;
        elseif path(count) == desired_point(count2)
            epsilon = 0.5; % How close zeta must be for the final switch, need more precision to get closer to the desired
            count2 = count2 + 1;
        else
            epsilon = 2; % How close zeta must be for intermediate switches
        end
        dwell_time = ceil(0.5*(log((abs(zeta-data(path(count),13))/epsilon)+1))/(log(1/delta_z))) % Dwell-time computation; see Prop 2, Veer et al CDC 2017
        k = path(count);
        last_step = i;
        count = count+1;
    end
    
    desired_speed = data(k,1);
    [betta,beta_c,theta_minus,theta_plus]=fcn_alpha_red_correction(ival',XX(:,i),N);
    %     [betta,beta_c,theta_minus,theta_plus]=fcn_alpha_red_correction(ival',fp_original,N);
    beta_ctrl = pinv(G)*data(k,2); % Compute the control parameters for k-th orbit
    beta_s_p = fcn_beta_correction(XX(:,1),XX(:,i),beta_ctrl);
    for j=1:1:4
        beta_s_p_tau(j,:) = fcn_coeffs_theta_to_tau(beta_s_p(j,:),theta_plus,theta_minus);
    end
    Fext = [0 ; 0];
    options = odeset('Events',@(t,q)touchdown5_steps(t,q,0),'RelTol',1e-5,'AbsTol',1e-4);
    [t,y,te,ye,ie] = ode45(@(t,Z) return_map5(t+tplot,Z,Fext,tF,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N),[0 20],XX(:,i),options); % Swing Phase
    simout.t=[simout.t;t+tplot];
    [pstk ph pswk pswf phead vh Menergy]=out_kinematics(y);
    
    simout.stk = [simout.stk  ; pstk];
    simout.h   = [simout.h    ; ph];
    simout.swk = [simout.swk  ; pswk];
    simout.swf = [simout.swf  ; pswf];
    simout.head= [simout.head ; phead];
    
    ntouch=ntouch+size(pswf,1);
    ptouch=ptouch+pswf(size(pswf,1),1);
    str_lng(i,1)=ntouch;
    str_lng(i,2)=ptouch;
    
    zeta = 0.5*(gamma0(ye)*ye(6:10)')^2; % Compute zeta
    
    XX_pi = ye';
    [XX(:,i+1) , Fimp] = impact_map5(ye); % Impact
    gam0_plus=gamma0(XX(:,i+1));
    gam0_minus=gamma0(ye);
    delta_zero=(gam0_plus*XX(6:10,i+1))/(gam0_minus*ye(6:10)'); % Compute delat_zero; see Theorem 1 of Veer et al CDC 2017
    timpact(i)=tplot;
    impactline=[timpact(i) timpact(i)];
    
%% Compute the input/outputs/friction cone for verification
    clear u h hdot hc hcdot h_original hdot_original GRF s friction_cone;
    clear output_error_d output_error_c output_error_s output_error
    for tp=1:size(t)
        [fx,gx,gforce]=fcn_stance_dynamics5(y(tp,:));
        u(tp,:)=fcn_stance_controller5_correction(y(tp,:)',fx,gx,gforce,Fext,betta,beta_c,beta_s_p_tau,theta_minus,theta_plus,M_bez,N);
        if M_bez==6
            [h(tp,:),hdot(tp,:),~]=fcn_liederivative(y(tp,:),betta,theta_minus,theta_plus,M_bez);
            [hc(tp,:),hcdot(tp,:),~]=fcn_liederivative_correction(y(tp,:)',beta_c,theta_minus,theta_plus,N,0.5);
            [hs(tp,:),hsdot(tp,:),~]=fcn_liederivative_correction(y(tp,:)',beta_s_p_tau,theta_minus,theta_plus,7,0.9);
            [h_original(tp,:),hdot_original(tp,:)]=fcn_liederivative(y(tp,:),betta,theta_minus,theta_plus,M_bez);
        elseif M_bez==5
            [h(tp,:),hdot(tp,:)]=fcn_liederivative_m5(y(tp,:),betta);
        end
        s(tp)=(c*y(tp,1:5)'-theta_plus)/(theta_minus-theta_plus);
        GRF(tp,:)=fcn_GRF(y(tp,:),u(tp,:));
        friction_cone(tp) = abs(GRF(tp,1)/GRF(tp,2));
        cm(tp,:)=fcn_cm_kinematics(y(tp,1:5));
        theta(tp) = c*y(tp,1:5)';
    end
    
    
    output_error_c=[hc hcdot];
    output_error_s=[hs hsdot];
    output_error_d=[h hdot];
    output_error=[(h-hc-hs) (hdot-hcdot-hsdot)];
    norm_output_error = max(max(abs(output_error)));
    norm_output_error_c = max(max(abs(output_error_c)));
    norm_output_error_s = max(max(abs(output_error_s)));
    norm_output_error_d = max(max(abs(output_error_d)));
    
    step_length(i,:)=fcn_position_swingfoot(ye);
    average_speed(i)=step_length(i,1)/te;
    
    if i>=2
        accel(i)=(average_speed(i)-average_speed(i-1))/te;
    end
    
    friction_cone = [friction_cone abs(Fimp(1)/Fimp(2))];
    fric_max = max(friction_cone);
    GRF_hor_min = min(GRF(:,1));
    [u_max,u_max_index] = max(max(abs(u')));
    u_max;
    % Print the desired speed, average speed, max input torque, friction cone, and maxium output error over that step
    fprintf('i=%1.0f,\t desired_speed=%1.4f,\t average_speed=%1.4f,\t umax= %3.2f,\t friction_cone= %2.4f \t max_output_error= %1.5f \n',i,desired_speed,average_speed(i),u_max,fric_max,max(max(output_error(:,:))));
    tplot=tplot+te;
    clear output_error norm_output_error hs hsdot theta;
    i = i+1;
end