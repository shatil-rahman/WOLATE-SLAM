function [ r_k_k_minus, C_k_k_minus, rho,exit_code] = WOLATE( r_jk, r_jk_minus)
%WOLATE Computes optimal C_k_kminus and r_k_k_minus 
    
    M = length(r_jk(1,:));
    b_k_alpha = zeros(3*M,1);
    b_k_beta = zeros(3*M,1);
    A_k_alpha = zeros(3*M, 6+3*M);
    A_k_beta = zeros(3*M, 6+3*M);
    W_inv = zeros(6*M,6*M);
    indexes = [];
    eps = 0.000001;
    
    r_k_k_minus.value = zeros(3,1);
    r_k_k_minus.cov = zeros(3,3);
    C_k_k_minus.value = eye(3);
    C_k_k_minus.cov = zeros(3,3); %The covariance is of the ROTATION VECTOR NOT THE DCM!!!!!
    rho_j_k(1,1:20) = struct('value',zeros(3,1),'cov',zeros(3,3));
    R_jk = r_jk.cov;
    R_jk_minus = r_jk_minus.cov;
    %% Create the A and B matrices
    for j = 1:M
        if isequal(r_jk_minus(:,j),zeros(3,1)) || isequal(r_jk(:,j),zeros(3,1))
            continue;
            
        end
        
        b_k_alpha(3*j-2:3*j) = -r_jk_minus(:,j);
        b_k_beta(3*j-2:3*j) = -r_jk(:,j);
        
        temp = [eye(3), cross_oper(r_jk(:,j) + r_jk_minus(:,j))];
        A_k_alpha(3*j-2:3*j, 1:6) = temp;
        A_k_alpha(3*j-2:3*j, 6 + 3*j-2:6+3*j) = -eye(3);
        A_k_beta(3*j-2:3*j, 6 + 3*j-2:6+3*j) = -eye(3);
%         W_j = R_jk + R_jk_minus - ((R_jk_minus - R_jk).'/(R_jk + R_jk_minus))*(R_jk_minus - R_jk);
        W_j = R_jk_minus - (R_jk_minus.'/(R_jk + R_jk_minus))*R_jk_minus;
        W_inv(3*j-2:3*j,3*j-2:3*j) = inv(W_j);
        W_inv(3*M + 3*j-2:3*M+ 3*j,3*M + 3*j-2:3*M + 3*j) = inv(R_jk);
        indexes = [indexes, 3*j-2: 3*j];
        
    end
    
    %%Check if any common landmarks
    
    if isempty(indexes)
        exit_code = -1;
        return
    end
        
    
    b_k_alpha = b_k_alpha(indexes);
    b_k_beta = b_k_beta(indexes);
    
    A_k_alpha = A_k_alpha(indexes,[1:6,6+indexes]);
    A_k_beta = A_k_beta(indexes, [1:6,6+indexes]);
    
    W_inv = W_inv([indexes, 3*M + indexes], [indexes, 3*M + indexes]);
    
    b_k = [b_k_alpha; b_k_beta];
    A_k = [A_k_alpha; A_k_beta];
  
   %%Least squares solution, add a tiny regularization term to make sure
   %%A.'W^-1A is spd
    
    x_k = ((A_k'*W_inv*A_k) + eps*eye(length(A_k(1,:))))\(A_k'*W_inv)*b_k;
    
    System_info = (A_k'*W_inv*A_k) + eps*eye(length(A_k(1,:)));
    
    %% Means
    
    r_k_k_minus.value = -(eye(3) + cross_oper(x_k(4:6)))\x_k(1:3);
    C_k_k_minus.value = (eye(3) + cross_oper(x_k(4:6)))\(eye(3) - cross_oper(x_k(4:6)));
    
    a = x_k(4:6)./(norm(x_k(4:6)));
    theta = 2*atan((norm(x_k(4:6))));
    rotation_vec = a.*theta;
    
    %% Covariances 
    
    
    
    r_k_k_minus.cov = (eye(3) + cross_oper(x_k(4:6)))/System_info(1:3,1:3)*(eye(3) + cross_oper(x_k(4:6))).';
    C_k_k_minus.cov = inv(System_info(4:6,4:6));
    rho = zeros(3*M,1);
    rho(indexes) = x_k(7:end);
    for j=1:M
       rho_j_k(1,j).value = rho(3*j-2:3*j);
       rho_j_k(1,j).cov = inv(System_info(6+3*j-2:6+3*j));
    end
%     rho_jk = x_k(7:end);
    exit_code = 1;

end

