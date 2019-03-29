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

function model = load_GSET(filename, rnd_seed, p_cut)

  if(nargin == 1)
    p_cut = -1;
  else
    %set problem generator seed (for x0)
    rng(rnd_seed);
  end
  

  delimiter = ' ';
  headLineN = 1;
  
  data_in = importdata(filename,delimiter,headLineN);
  n = str2num(data_in.textdata{1});
  N = n(1);
  n_nnz = n(2);
  Q = zeros(N,N);
  for ii = 1: n_nnz
    row = data_in.data(ii,1);
    col = data_in.data(ii,2);
    val = data_in.data(ii,3); %maximization problem!
    if(row ~= col) %should not be the case, but...
      Q(row,col) = val;
      Q(col,row) = val;
    end
  end
  if(nnz(Q) ~= n_nnz*2)
    error('Error. Something is wrong with the input file, nnz wrong')
  end
  d = sum(Q,2);
  %Q = Q - (diag(d) + 1e-20);
  Q = Q - diag(d);
  
  x0 = zeros(N,1);
  if(p_cut > 0)
    k_cut = p_cut * N;
    if(k_cut > N)
      error('Error. Cut larger than size');
    end
    Aeq = ones(1,N);
    beq = k_cut;
    d=datasample([1:N], k_cut,'Replace',false);  
    x0(d) = 1;
  else    
    Aeq = sparse(0,N);
    beq = [];
  end
  binary=[1:N]';
  lb = zeros(N,1);
  ub = ones(N,1);
  Aeq_g = Aeq;
  beq_g = beq;
  Aeq_l = sparse(0,N);
  beq_l = [];

  model.size = N;
  model.Q = sparse(Q);
  model.c = sparse(N,1);
  model.Aeq = sparse(Aeq_g);
  model.beq = beq_g;  
  model.Aineq = sparse(0,N);
  model.bineq = [];
  model.x0 = x0;

  model.integers = [];
  model.binary = binary;

  model.lb = lb;  
  model.ub = ub;
  model.const = 0;
  model.local_constraints.Aeq = sparse(Aeq_l);
  model.local_constraints.beq = beq_l; 
  model.local_constraints.Aineq = sparse(0,N);
  model.local_constraints.bineq = [];
  model.local_constraints.lb=lb;  
  model.local_constraints.ub=ub;


end
