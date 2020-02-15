function sol = rp_admm_lasso(X, y, block_num, alpha, lambda_penalty, lambda_admm, max_iter)

%  RP ADMM without prefactorization
%
%  Solver for following Elastic Net Problem 
%  
%  min (1/2N)* (y-X*beta)^T(y-X*beta) + lambda( (1-alpha)/2*beta'beta + (1-alpha)*||beta||_1 )
%
%  alpha elastic net parameter; lambda penalty of elastic net; gamma penalty of ADMM
%
%  Written for MATLAB_R2016b
%  Copyright (C) 2019
%     Kresimir Mihic <kmihic@alumni.stanford.edu>
%     Mingxi Zhu <mingxiz@stanford.edu>
%     Yinyu Ye <yyye@stanford.edu>

rac_time_start = tic;
rac_preparing_time_start = tic;

[num_obs, num_features] = size(X);
x0 = zeros(num_features,1);
x_k = x0;
x_k_2 = x0;
y_k = zeros(num_features,1);

% prepare and scale
c = -1/num_obs*X'*y;
A = 1/sqrt(num_obs)*X;

tol_temp = inf;
sub_model_time = 0;
sub_model_time_2 = 0;
sub_model_time_3 = 0;
sub_model_time_4 = 0;
sub_solver_time = 0;
num_iter = 0;

% x0 = 0, Ax = 0;
Ax = zeros(num_obs, 1);


%RP prepare sub matrix
order = randperm(num_features);
block_size_r = floor(num_features/block_num);
block_size_n = num_features-(block_num - 1)*block_size_r;
store_A_Sub_A_sub=cell(1,block_num);
for i_block = 1:1:block_num
    if i_block == block_num
        A_sub =  A(:,order(((block_num - 1)*block_size_r + 1):num_features));
        store_A_Sub_A_sub{i_block} = (A_sub'*A_sub);
    else
        A_sub = A(:,order(((i_block - 1)*block_size_r + 1):i_block*block_size_r));
        store_A_Sub_A_sub{i_block} = (A_sub'*A_sub);
    end
end



rac_preparing_time = toc(rac_preparing_time_start);



while num_iter<max_iter

%submodeltime for c_dual_vector
sub_model_time_2_start = tic;

c_res = c - y_k;
x_block_index_perm = randperm(block_num);
%x_index_perm = 1:dim;
x_visited_block = 1;

sub_model_time_2 = toc(sub_model_time_2_start) + sub_model_time_2;
%update primal
  while x_visited_block <= block_num
  
  sub_model_time_start = tic;
  %random select block num from x_unvisited [1,num_features] 
  current_block_index = x_block_index_perm(x_visited_block);
  if current_block_index == block_num
          x_update_index = order(((block_num - 1)*block_size_r + 1):num_features);
          block_size = block_size_n;
      else
          x_update_index =  order(((current_block_index - 1)*block_size_r + 1):current_block_index*block_size_r);
          block_size = block_size_r;
  end
  
  
  %model buildup
  A_sub = A(:,x_update_index);
  H_current = store_A_Sub_A_sub{current_block_index};
  c_current_sub = A_sub'*Ax - H_current*x_k(x_update_index);
  H_current = H_current + lambda_admm(num_iter+1)*speye(block_size);
  c_current = c_current_sub + c_res(x_update_index) - lambda_admm(num_iter+1)*x_k_2(x_update_index); 
  
  
  %model solver LDL factorization
  %% min x^TQx+c^Tx
  %% FOC (2Q+betaA^TA)x=-c
  right_side = -c_current;
  left_side_unfactorized = H_current;
  sub_model_time = toc(sub_model_time_start) + sub_model_time;
  
  sub_solver_time_start = tic;
  results_x = left_side_unfactorized\right_side; 
  sub_solver_time = toc(sub_solver_time_start) + sub_solver_time;
  
  %update Ax
  sub_model_time_3_start = tic;
  diff_x = x_k(x_update_index) - results_x;
  Ax = Ax - A_sub*diff_x;
  
  %update primal
  x_k(x_update_index) = results_x;
    
  %update x_k_2
  penalty_co = 1/((1-alpha)*lambda_penalty(num_iter+1)+lambda_admm(num_iter+1));
  x_k_2(x_update_index) = penalty_co*soft( y_k(x_update_index) - lambda_admm(num_iter+1)*x_k(x_update_index) ,lambda_penalty(num_iter+1)*alpha);   
  
  
  
  x_visited_block = x_visited_block + 1;
  sub_model_time_3 = toc(sub_model_time_3_start) + sub_model_time_3;

  
  end
  
%update dual
  sub_model_time_4_start = tic;  
  res_k = (x_k - x_k_2);
  
%update dual y_k,z_k
  y_k = y_k - lambda_admm(num_iter+1)*res_k;
  
  num_iter = num_iter+1;
    
  sub_model_time_4 = toc(sub_model_time_4_start) + sub_model_time_4;
  

end 

rac_time = toc(rac_time_start);

rac_admm_x = x_k;
sub_model_time_f = sub_model_time + sub_model_time_2 + sub_model_time_3 + sub_model_time_4;

sol.beta = rac_admm_x;
sol.beta2 = x_k_2;
sol.total_time = rac_time;
sol.solver_time = sub_solver_time;
sol.model_time = sub_model_time_f;
sol.model_detail = [sub_model_time, sub_model_time_2, sub_model_time_3, sub_model_time_4];
sol.prepare_time = rac_preparing_time;
sol.num_iter = num_iter;
sol.prim_tol = tol_temp;