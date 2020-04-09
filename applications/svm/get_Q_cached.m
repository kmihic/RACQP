 %
 % RACQP -  Randomly Assembled Cyclic ADMM Quadratic Programming Solver
 % Copyright (C) 2019 
 %     Kresimir Mihic <kmihic@alumni.stanford.edu>
 %     Mingxi Zhu <mingxiz@stanford.edu>
 %     Yinyu Ye <yyye@stanford.edu>
 %
 % This file is part of RACQP 
 %
 %
 
function Q = get_Q_cached(x_i,x_j,x1,x2)
 
  global X_train;
  global X_train_sq;
  global Y_train;
  global Q_cached;
  global is_cached;
  global N_x;
  global M_x;
  global Col_cnt;
 
  % hard coding to use 70% free memory
  P_CACHE = 0.7; 
  USE_SWAP = false;

  if(isempty(is_cached))
    is_cached = zeros(N_x,1);
    %check the max Q size we can fit in
    % The matrix is double, thus 8 bytes per value
    % So, the total size is n*m*8+overhead
    free_mem = get_free_memory(USE_SWAP);
    M_x = min(N_x,floor(free_mem*P_CACHE/(N_x*8)));
    Q_cached = zeros(N_x,M_x);
    Col_cnt = 0;
    disp("Using cached Q. Mem: "+(1e-9*free_mem)+" Gb. Max number of columns: "+M_x)
  end
 
  % x'Qx used only for objVal cals. Not needed for SVM
   if(nargin == 4)
     Q = NaN;
     return;
   end
 
  % in the code x_i is usually Q(':',x_j), except for when RP mode is used
  % and then Q(x_ix,x_ix) is called only once. No reason to complicate the code
  % if no (':',x_j), just calculate the kernel and get back 
  if(strcmpi(x_i,':'))
    % Q(:,:) is called only twice, from get_residuals and obj_val at the end of rac
    % returning an empty array, and skipping residual calc
    if(strcmpi(x_j,':'))
      Q = sparse(N_x,N_x);
      return;
    end
    x_ix = is_cached(x_j);
    r_ix = x_j(x_ix==0);   
    % construct missing
    if(length(r_ix > 0))
      K = bsxfun(@minus,X_train_sq(x_i),(2*X_train(x_i,:))*X_train(r_ix,:).');
      K = bsxfun(@plus,X_train_sq(r_ix).',K);
      K = exp(-K);      
      q = bsxfun(@times, Y_train(x_i), K);
      q = bsxfun(@times, q, Y_train(r_ix)');
      % save the columns if space left
      if(Col_cnt < M_x)
        be = min(M_x,Col_cnt+length(r_ix));
        ix = [(Col_cnt+1):be];
        j_ix = 1:length(ix);
        %Q_cached(:,ix) = exp(-K_r(:,j_ix));
        Q_cached(:,ix) = q(:,j_ix);
        saved_cols = x_j(j_ix);
        is_cached(saved_cols) = ix;
        Col_cnt = be;
      end
      p_rix = find(ismember(x_j,r_ix));
      Q(:,p_rix) = q;
    end
    
    % add stored
    x_cached = x_ix>0;
    c_ix = x_ix(x_cached);
    if(length(c_ix > 0))
      q_ix = Q_cached(:,c_ix);
      Q(:,x_cached) = q_ix;
    end
  else 
    K = bsxfun(@minus,X_train_sq(x_i),(2*X_train(x_i,:))*X_train(x_j,:).');
    K = bsxfun(@plus,X_train_sq(x_j).',K);
    K = exp(-K);
    Q = bsxfun(@times, Y_train(x_i), K);
    Q = bsxfun(@times, Q, Y_train(x_j)');
  end
 
  
  if(nargin == 3)
    Q = Q*x1;
  end
 
%end
