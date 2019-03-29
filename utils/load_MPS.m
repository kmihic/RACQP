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


function model=load_MPS(filename, diag_Q, sparsity_Q)
  
  s=mpsread(filename);
  model.size = length(s.f);
  if(diag_Q)
    model.Q = speye(model.size,model.size);
  elseif(nargin > 2 )
    p_setup.n_var = model.size;
    p_setup.Q_sparsity = sparsity_Q;
    model.Q = getQ(p_setup);
  else
    model.Q = sparse(model.size,model.size);
  end
  model.c = sparse(s.f);
  model.Aeq = sparse(s.Aeq);
  model.beq = s.beq;
  model.Aineq = sparse(s.Aineq);
  model.bineq = s.bineq;
  model.x0 = zeros(model.size,1);
  model.integers = s.intcon;
  model.binary = [];
  model.lb = s.lb;
  model.ub = s.ub;
  model.const = 0;
end

function Q=getQ(n_var, sparsity)
  if(sparsity==1)
    Q=sparse(n_var,n_var);
    return;
  end
  A=sprand(n_var,n_var,1-sparsity);
  Q = A+A';
  %Q has all positive values (rand>=0), no need for abs
  %make it psd
  off_diag = sum(Q,2) - diag(Q) +1e-10; 
  max_off_diag = max(off_diag);
  Q = Q+diag(off_diag);
  %normalize Q
  Q = Q./max_off_diag;
  Q = sparse(Q);
end