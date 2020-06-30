function beta_s_p_tau = fcn_hs_theta_to_tau(beta_s_p,theta_plus,theta_minus)
%FCN_HS_THETA_TO_TAU Summary of this function goes here
%   Detailed explanation goes here
for m=1:1:4
    beta_s_p_tau(m,:) = fcn_coeffs_theta_to_tau(beta_s_p(m,:),theta_plus,theta_minus);
end

end

