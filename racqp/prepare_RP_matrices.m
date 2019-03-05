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

function matrices = prepare_RP_matrices(Q, A,groups, beta, factorize, ...
                    use_sparse, do_xk)

  matrices={};
  n_blocks = size(groups,1);
  beta_h = beta/2;
  for block_ix = 1:n_blocks
    x_ix = groups{block_ix,1};
    A_sub = A(:,x_ix);
    Q_current = Q(x_ix,x_ix)+beta_h*A_sub'*A_sub;
    if(do_xk)
      Q_fact = Q_current+beta_h*speye(length(x_ix));
    else
      Q_fact = Q_current;
    end
    if(~factorize)
      R=[];
      S=[];
    elseif(use_sparse)
        [R, p, S]=chol(2*Q_fact);
    else
      R = chol(2*Q_fact);
      S = [];
    end
    matrices(block_ix,:) = {Q_current, R, S, x_ix};
  end


end

