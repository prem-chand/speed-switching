function [betta,theta_minus,theta_plus,event]=fcn_transition_ctr(fixp1,fixp2)

[betta1,theta_minus_1,theta_plus_1]=fcn_alpha_red(fixp1);
[betta2,theta_minus_2,theta_plus_2]=fcn_alpha_red(fixp2);
M1=(max(size(fixp1))-10)/4+2;
M2=(max(size(fixp2))-10)/4+2;
M12=min(M1,M2);

event=fixp2(2);
theta_minus = theta_minus_2;
theta_plus  = theta_plus_1;
betta(:,1)= betta1(:,1);
betta(:,2)= betta1(:,1) - M1/M12*(theta_minus_2 - theta_plus_1)/(theta_minus_1 - theta_plus_1)*(betta1(:,1)-betta1(:,2));
%betta(:,3:M-1)=(betta1(:,3:M-1)+betta2(:,3:M-1))/2;
%betta(:,3:M-1)=(betta1(:,3:end-2)+betta2(:,3:end-2))/2;
if M1==M2
    betta(:,3:M12-1)=0.5*(betta1(:,3:M1-1)+betta2(:,3:M2-1));
elseif M1==M12
    betta(:,3:M12-1)=betta1(:,3:M1-1);
elseif M2==M12
    betta(:,3:M12-1)=betta2(:,3:M2-1);
end

%betta(:,3:M)=(betta1(:,3:M)+betta2(:,3:M))/2;
betta(:,M12)= betta2(:,end) + M2/M12*(theta_minus_2 - theta_plus_1)/(theta_minus_2 - theta_plus_2)*(betta2(:,end-1)-betta2(:,end));
betta(:,M12+1)= betta2(:,end);
end