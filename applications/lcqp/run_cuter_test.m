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

% Verify results: CUTEr benchmark problems
%      Loading data from OSQP solver dataset (mat files)
%      Original mps files:
%      http://www.cuter.rl.ac.uk/Problems/marmes.html
%

function run_cuter_test(solver,r_time,epsilon, max_iter, rnd_seed)

  inst = get_instances();
  solutions = [];
  data_path = "../data/data_cuter/";
  for ii = 1:size(inst,1)
    filename = data_path+inst(ii,1);
    disp('Solving: '+filename);
    disp("LOADING THE MODEL...")
    load(filename);
    % Mosek complains, other solvers not, thus to make all run smoothly
    % without building a cubic model for Mosek (as per suggestions)
    % check eigvals if positive, and if not adjust
    % by adding a small diagonal.
    if( str2num(inst(ii,4))~= 0 && strcmpi(solver,'mosek'))
       scale_Q = norm(model.Q,Inf);
       model.Q = model.Q/scale_Q;
       model.c = model.c/scale_Q;
       model.Q = model.Q + speye(model.size) * str2num(inst(ii,4));
    end
    if(strcmpi(solver,'racqp'))
       m_orig = model;
       model.x0 = get_init_pt(model.lb, model.ub);
       [model,max_Q,max_Aeq, max_Aineq] = scale_model(model, true,true);
       model.scale_Aeq = max_Aeq;
       model.scale_Aineq = max_Aineq;
       model.norm_coeff = max_Q;
       %clean the model from embedded, unbounded slacks
       model = clean_model_constraints(model, true);
    end
     inst_param = get_instance_run_params(model, solver, r_time, false, ...
                          epsilon, max_iter);
     inst_param.racqp_run_p = get_cuteR_run_params(str2num(inst(ii,2)), inst(ii,3), ...
                     rnd_seed,epsilon, max_iter, r_time);
     s = solve_instance(inst_param);
     s.name = inst(ii,1);
     solutions = [solutions,s];
  end
  print_solutions(solutions);
end

function run_params = get_cuteR_run_params(beta, pL_mode, rnd_seed,epsilon,max_iter,max_time)

  %load default runtime parameters
  n_blocks = 1;
  run_params = default_run_params(n_blocks, beta,max_time,epsilon,max_iter);
  %change some parameters
  run_params.rnd_seed = rnd_seed;
  %turn on verbose
  %run_params.debug = 1;

  if(strcmpi('l_eq', pL_mode))
    run_params.single_block_solver = 'diag_kkt';
  elseif(strcmpi('l_bounds',pL_mode))
    run_params.single_block_solver = 'bounds_by_gurobi';  
    gurobi_params = default_gurobi_parameters;
    %change max time for gurobi for each subproblem
    gurobi_params.TimeLimit = 200;     
    run_params.gurobi_params = gurobi_params;
  elseif(strcmpi('none', pL_mode))
    run_params.single_block_solver = 'cholesky';
  elseif(strcmpi('multi', pL_mode))
    run_params.n_blocks = 2;
  end
end

function inst = get_instances()

inst=[
  "AUG2DC","1","l_eq","0";
  "AUG2DQP","50","l_eq", "0";
  "AUG3DC","1","l_eq", "0";
  "AUG3DCQP","1","l_eq", "0";
  "BOYD1","1","l_eq", "0";
  "CONT-050","350","l_eq", "0";
  "CONT-100","350","l_eq", "0";
  "CONT-101","350","l_eq", "0";
 "CONT-300","350","l_eq", "0";
 "CVXQP1_L","1","none","1e-16";
  "CVXQP2_L","1","none","1e-16";
  "DTOC3","25","l_eq", "0";
  "HUES-MOD","1","l_eq", "0";
  "HUESTIS","1","l_eq", "0";
  "UBH1","12000","l_bounds", "0"];


end


