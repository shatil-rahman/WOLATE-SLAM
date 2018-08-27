%% Script to generate data for a bundle adjustment problem
clear all
close all
clc;

%%
load dataset3.mat

k_start = 970;
k_end = 1100;
timesteps = k_start:k_end;
noise = 0.01;

measurement_lists = gen_meas(timesteps, noise);
K = length(measurement_lists(:,1));

r_ka_chained = zeros(3,K);
r_ja_pseudo = zeros(3,K,20);
r_ja_pseudo(1:K,1:20) = struct('meas',zeros(3,1),'cov',zeros(3,3));
% r_ka_chained(:,1) = r_i_vk_i(:,timesteps(1));
% C_k_minus_a = ax_ang2dcm(theta_vk_i(:,timesteps(1)));
C_k_minus_a = eye(3);
r_ja_pseudo(1,:) = measurement_lists(1,:);
% r_pseudo_cov(:,1,:) = [noise;noise;noise];


for k = 2:K
    meas_k = measurement_lists(k,:);
    meas_k_minus = measurement_lists(k-1,:);
    [r_k_kminus, C_Kk, rho_est,flag] = WOLATE(meas_k,meas_k_minus);
    
    if isempty( r_k_kminus)
        error('No common landmarks between pose k-1 and k')
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

[ r_ka_a,C_ka,r_j_a_est] = BA_WOLATE( r_ja_pseudo(:,1:end,:), measurement_lists(:,1:end,:));
r_ka_a_true = zeros(3,K);
r_ja_a_true = zeros(3,20);

%%
for k = 1:K
   r_ka_a_true(:,k) = r_i_vk_i(:,timesteps(k)) - r_i_vk_i(:,k_start);
   r_ka_a_true(:,k) = ax_ang2dcm(theta_vk_i(:,k_start))*r_ka_a_true(:,k);
   
end

for j = 1:20
   r_ja_a_true(:,j) = rho_i_pj_i(:,j) - r_i_vk_i(:,k_start);
   r_ja_a_true(:,j) = ax_ang2dcm(theta_vk_i(:,k_start))*r_ja_a_true(:,j);    
end

%% Aligning the results

% [params,r_ka_a_aligned] = absor(-r_ka_a, r_ka_a_true(:,2:end));
[params2,r_j_a_est_aligned] = absor(r_j_a_est, r_ja_a_true);
for k = 1:K-1
   r_ka_a_aligned(:,k) = params2.R*(r_ka_a(:,k)) + params2.t;
end
%% Plotting

plot3(r_ka_a_true(1,1:end),r_ka_a_true(2,1:end),r_ka_a_true(3,1:end))
hold on
grid on
plot3(r_ka_a(1,:),r_ka_a(2,:),r_ka_a(3,:))
scatter3(r_ja_a_true(1,:), r_ja_a_true(2,:), r_ja_a_true(3,:))
scatter3(r_j_a_est(1,:), r_j_a_est(2,:), r_j_a_est(3,:))
plot3(r_ka_chained(1,2:end),r_ka_chained(2,2:end),r_ka_chained(3,2:end))
legend('true','BA','true','BA','daisy-chained WOLATE')

%% Error Calculations




