%% Script to generate data for a bundle adjustment problem
clear all
close all
clc;

%%
load dataset3.mat

timesteps = 1000:1100;

measurement_lists = gen_meas(timesteps, 0.00001);
K = length(measurement_lists(1,:,1));

r_ka_chained = zeros(3,K);
r_ja_pseudo = zeros(3,K,20);
% r_ka_chained(:,1) = r_i_vk_i(:,timesteps(1));
% C_k_minus_a = ax_ang2dcm(theta_vk_i(:,timesteps(1)));
C_k_minus_a = eye(3);
r_ja_pseudo(:,1,:) = gen_meas(timesteps(1),0.00001);


for k = 2:K
    meas_k = squeeze(measurement_lists(:,k,:));
    meas_k_minus = squeeze(measurement_lists(:,k-1,:));
    [r_k_kminus, C_Kk, rho_est,flag] = WOLATE(meas_k,meas_k_minus, eye(3).*0.00001,eye(3).*0.00001);
    
    if isempty( r_k_kminus) 
        break
    end
    C_k_a = C_Kk*C_k_minus_a;
    r_ka_chained(:,k) = r_ka_chained(:,k-1) + C_k_a.'*r_k_kminus;
    for j = 1:20
        if isequal(rho_est(:,j),zeros(3,1))
            continue;
        end
        r_ja_pseudo(:,k,j) = r_ka_chained(:,k) + C_k_a.'*rho_est(:,j);
    end
    C_k_minus_a = C_k_a;
    
end

%% Bundle Adjustment WOLATE

[ r_ka_a,C_ka,r_j_a_est] = BA_WOLATE( r_ja_pseudo(:,2:end,:), gen_meas(timesteps(2:end),0.00001));


