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

%
%TEST: Markowitz portfolio selection and min variance problem
%      Data source Center for Research in Security Price (CRSP)
%      Time period Jan-Dec 2018
%
%

function [solutions inst_param] = run_portfolio_test(solver, r_time, low_rank,binary, ...
                  epsilon, max_iter, rnd_seed, quiet, inst)

  if(nargin <= 4)
    epsilon = 1e-5;
  end
  if(nargin <= 5)
    max_iter = 4000;
  end
  if(nargin <= 6)
    rnd_seed = 123;
  end
  if(nargin <= 7)
    quiet = false;
  end
  if(nargin <= 8)
    if(low_rank)
      inst = get_instances_lr();
    else    
      inst = get_instances_reg();
    end
  end
   
  %setup problem data
  kappa = 1e-5;
  if(binary)
    p_cut = 0.5;
  else
    p_cut = -1;
  end
  solutions = [];
  data_path = "../data/data_markowitz/";
  for ii = 1:size(inst,1)
    filename = data_path+inst(ii,1);
    disp('Solving: '+filename);
    disp("LOADING THE MODEL...")
    model = load_Markowitz_model(filename, p_cut, kappa, low_rank, binary);
    inst_param = get_instance_run_params(model, solver, r_time, binary, ...
                         epsilon, max_iter);
    if(strcmpi(solver,'racqp'))
      if(binary)
        inst_param.racqp_run_p = get_portfolio_run_params(model.size, r_time,low_rank);
      else
        inst_param.racqp_run_p = get_min_variance_run_params(model.size, r_time,low_rank,...
                      epsilon, max_iter, rnd_seed);
      end
    end
    s = solve_instance(inst_param);
    s.name = inst(ii,1);
    % opt/best obj_val is for binary problems
    if(binary)
      s.obj_val = inst(ii,2);
    end
    solutions = [solutions,s];
  end
  if(~quiet)
    if(binary)
      print_solutions_binary(solutions);
    else
      print_solutions(solutions);
    end
  end
end

function run_mip = get_portfolio_run_params(nsize, max_time, low_rank)
  
  if(low_rank)
    n_blocks = 50;
    beta = 0.5;
  else
    n_blocks = 100;
    beta = 0.05;
  end
  lambda=0.4;
  p_trial = 0.005;
 
  %load default runtime parameters
  run_params = default_run_params(n_blocks, beta,max_time);
  %change some parameters
  %turn on verbose
  %run_params.debug = 1;
 
  %load default gurobi parameters
  gurobi_params = default_gurobi_parameters;
  %change how long gurobi will spend on each subproblem
  gurobi_params.TimeLimit = 2;
  run_params.gurobi_params = gurobi_params;
  
  %set MIP run parameters
  run_mip.permute_type = 'swap';
  run_mip.permute_dist = 'exponential'; 
  run_mip.permute_min = 2; 
  run_mip.permute_max = nsize; 
  run_mip.permute_mu = lambda*nsize; 
  run_mip.max_iter = Inf;
  run_mip.max_nperturb = Inf;
  run_mip.max_rtime = max_time; 
  run_mip.rnd_seed = 123; 
  run_mip.n_perturb_trial = max(2,round(p_trial*nsize));
  run_mip.run_sub = run_params;
  run_mip.debug = 0;
  run_mip.mip_epsilon = 1e-5;
  
end

  
function run_params = get_min_variance_run_params(nsize, max_time, low_rank, ...
                     epsilon, max_iter, rnd_seed)
  
  %setup problem data
  if(low_rank)
    n_blocks = 1;
  else
    n_blocks = 50;
  end
  beta = 1;
 
  %load default runtime parameters
  run_params = default_run_params(n_blocks, beta,max_time,epsilon,max_iter);
  %change some parameters
  run_params.rnd_seed = rnd_seed;
  %turn on verbose
  run_params.debug = 0;
 
  if(low_rank)
    %if using lowrank model, we use single block, kkt mode
    %with localized constraints. Q is diagonal, so model is sparse
    run_params.n_blocks = n_blocks;
    run_params.single_block_solver = 'diag_kkt';
    run_params.use_sparse = true;
  else
    %regular model has dense Q
    run_params.use_sparse = false;
  end
  run_params.n_blocks = n_blocks;
end
  
  
function model = load_Markowitz_model(filename, p_cut, kappa, low_rank, mip)
  
  load(filename);
  %get the problem 
  binary=[];
  if(low_rank)  
    ic=length(index_continuous);
    c = c_lowrank;
    N=length(c); 
    Nx=N-ic;
    Q = [kappa*speye(Nx),sparse(Nx,ic);sparse(ic,Nx),speye(ic)];
    beq = b_lowrank;
    Aeq = A_lowrank;  
    %x vars are either binary or boxed cont
    lb = zeros(N,1);
    lb(index_continuous) = -inf;    
    ub = ones(N,1);
    ub(index_continuous) = inf;  
    if(mip)
      n_x = length(index_integer);
      beq(1) = min(n_x,ceil(p_cut*n_x));    
      binary = index_integer;        
    else
      beq(1) = 1;
    end
  else
    N=size(Q,1);
    Q=sparse(Q) + kappa*speye(N);
    %x vars are either binary or boxed cont
    lb = zeros(N,1);  
    ub = ones(N,1);
    Aeq = sparse(ones(1,N)); 
    if(mip)
      beq = min(N,ceil(p_cut*N));
      binary = (1:N)';    
    else
      beq = 1;    
    end
  end
  
  
  model.size = N;
  model.c = c;
  model.Q = Q;
  model.Aeq = Aeq;
  model.beq = beq;  
  model.Aineq = sparse(0,N);
  model.bineq = [];
  model.x0 = zeros(N,1);
a = -100;
b = 100;
model.x0 = (b-a).*rand(model.size,1) + a;
  model.integers = [];
  model.binary = binary;
  model.lb = lb;  
  model.ub = ub;
  model.const = 0;
  if(mip)
    model.local_constraints.lb=model.lb;
    model.local_constraints.ub=model.ub;
    model.local_constraints.Aineq = sparse(0,N);
    model.local_constraints.bineq = [];
    model.local_constraints.Aeq = sparse(0,N);
    model.local_constraints.beq = [];
  end
end

function inst = get_instances_reg()
inst=[
"regular_quarterly","0.055";
"regular_monthly","0.144";
"regular_daily","1.164";  
];
end

function inst = get_instances_lr()
inst=[
"lowrank_quarterly","0.015"
"lowrank_monthly","0.104";  
"lowrank_daily","1.140";  
];
end
