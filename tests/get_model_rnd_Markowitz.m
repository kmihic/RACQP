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


function model = get_rnd_Markowitz(size, sparsity, mu, q_params, rnd_seed, num_chosen)

  %QP params
  p_setup.Q_sparsity = sparsity;
  p_setup.n_var = size;
  p_setup.c_sparsity = sparsity;
  p_setup.Aeq_sparsity = 1;
  p_setup.Aineq_sparsity = 1;
  p_setup.Aineq_n_row = 0;
  p_setup.Aeq_n_row = 0;
  p_setup.rnd_seed = rnd_seed;

  model = get_model_rnd_QP(p_setup,false,q_params);
  model.Q = model.Q + mu*speye(size);
  model.Aeq = sparse(ones(1,size)); 
  model.lb = zeros(size,1);  
  model.ub = ones(size,1);

  if(nargin > 5 && num_chosen > 0)
    if(num_chosen > size)
      error('Error. Cardinality must be [1,size)')
    end
    model.beq = num_chosen;
    model.binary = (1:size);   
  else
    model.beq = 1;  
  end
