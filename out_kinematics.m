function [pstk,ph,pswk,pswf,phead,vh,Menergy] = out_kinematics(x)

for i = 1:size(x,1)
        
        pstk(i,:)  = fcn_position_stanceknee(x(i,:)); 
        ph(i,:)    = fcn_position_hip(x(i,:));
        pswk(i,:)  = fcn_position_swingknee(x(i,:));
        pswf(i,:)  = fcn_position_swingfoot(x(i,:));
        phead(i,:) = fcn_position_head(x(i,:));
        vh(i,:)    = fcn_velocity_hip(x(i,:));
        
        Menergy(i,:)= fcn_energy(x(i,:));
end
end