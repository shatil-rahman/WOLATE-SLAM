function [ r_ka_a,C_ka,r_j_a_est] = BA_WOLATE( r_ja_k, r_jk_k )
%BA_WOLATE Summary of this function goes here
%   Detailed explanation goes here

    K = length(r_ja_k(1,:,1));
    M = 20;

    B = zeros(nnz(r_ja_k),1);
    A = zeros(length(B),6*K + 3*M);
    W_inv = zeros(length(B),length(B));
    
    
    meas_no = 1;
    for k = 1:K
       for j = 1:M
          
           if isequal(r_ja_k(:,k,j),zeros(3,1))
              continue 
           else
              B(meas_no:meas_no+2) = -r_jk_k(:,k,j);
              A(meas_no:meas_no+2, 6*k-5:6*k) = [eye(3), cross_oper(r_ja_k(:,k,j) +r_jk_k(:,k,j))];
              A(meas_no:meas_no+2, 6*K + 3*j-2: 6*K + 3*j) = -eye(3);
%               W_j = R_jk_minus - (R_jk_minus.'/(R_jk + R_jk_minus))*R_jk_minus;
%               W_j = r_jk_k(:,k,j) - (r_jk_k(:,k,j).'/(r_ja_k(:,k,j) + r_jk_k(:,k,j)))*r_jk_k(:,k,j);
%               W_inv(meas_no:meas_no+2,meas_no:meas_no+2) = eye(3);
              W_inv(meas_no:meas_no+2,meas_no:meas_no+2) = eye(3).*inv(k/nnz(r_ja_k(:,k,j)));
              
%               meas_no = meas_no + 3;
%               B(meas_no:meas_no+2) = -r_ja_k(:,k,j);
%               A(meas_no:meas_no+2, 6*K + 3*j-2: 6*K + 3*j) = -eye(3);
              
              meas_no = meas_no + 3;
              
           end
           
       end
    end
    System_mat = A.'*W_inv*A;
    System_mat = System_mat(7:end,7:end);
    Sys_vec = A.'*W_inv*B;
    Sys_vec = Sys_vec(7:end);
    
%     x = (A.'*A)\A.'*B;
    x  =System_mat\Sys_vec;
%     x = [0;0;0;x];
    
    r_ka_a = zeros(3,K-1);
    C_ka = zeros(3,3,K-1);
    r_j_a_est = reshape(x(6*(K-1)+1:end),3,M);
    
    
    
    for k = 1:K-1
        r_ka_a(:,k) = (eye(3) + cross_oper(x(6*k-2:6*k)))\x(6*k-5:6*k-3);
        C_ka(:,:,k) = (eye(3) + cross_oper(x(6*k-2:6*k)))\(eye(3) - cross_oper(x(6*k-2:6*k)));
    end

end

