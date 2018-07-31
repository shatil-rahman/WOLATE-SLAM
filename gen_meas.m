function [ rk_j_k, noise ] = gen_meas(timesteps, noise)
%GEN_MEAS Generates measurements in the vehicle frame
%   Detailed explanation goes here
    load dataset3.mat;
    rk_j_k = zeros(3,length(timesteps),20);
    for i=1:length(timesteps)
       k = timesteps(i);
       for j=1:20
          if(y_k_j(:,k,j)) ~= -1*ones(4,1)
              meas_inertial = rho_i_pj_i(:,j) - r_i_vk_i(:,k);
              phi_k = theta_vk_i(:,k);
              C_vk_i = ax_ang2dcm(phi_k);
              meas_jk = C_vk_i*meas_inertial;
              rk_j_k(:,i,j) = meas_jk + noise*randn(3,1);
          end
       end
    end

end

