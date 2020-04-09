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

function Q = get_Q(x_i,x_j,x1,x2)

 global Y_train;
 global X_train;
 global X_train_sq;
 %global g_gamma;

 % x'Qx used only for objVal cals. Not needed for SVM
  if(nargin == 4)
    Q = NaN;
    return;
  end
   % Q(:,:) is called only twice, from get_residuals and obj_val at the end of rac
   % returning an empty array, and skipping residual calc
  if(strcmpi(x_j,':'))
      Q = sparse(N_x,N_x);
      return;
   end

  K = bsxfun(@minus,X_train_sq(x_i),(2*X_train(x_i,:))*X_train(x_j,:).');
  K = bsxfun(@plus,X_train_sq(x_j).',K);
  K = exp(-K);
  Q = bsxfun(@times, Y_train(x_i), K);
  Q = bsxfun(@times, Q, Y_train(x_j)');

 if(nargin == 3)
   Q = Q*x1;
 end

end
