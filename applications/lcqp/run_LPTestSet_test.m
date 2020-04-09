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

% Verify results: LPTestSet benchmark problems + diag Q
%      http://plato.asu.edu/ftp/lptestset/


function solutions = run_LPTestSet_test(solver,r_time, epsilon, ...
              max_iter, rnd_seed,quiet, inst)

  if(nargin <= 5)
    quiet = false;
  end
  if(nargin <= 6)
    inst = get_instances();
  end
  solutions = [];
  data_path = "../data/data_lpTestSet/";
  for ii = 1:size(inst,1)
    filename = inst(ii,1);
    disp('Solving: '+data_path+filename);
    disp("LOADING THE MODEL...")
    load(data_path+filename);
    %add local bounds
    if(strcmpi(inst(ii,4),"l_bounds"))
      local_constraints.Aeq = sparse(0,model.size);
      local_constraints.beq = [];
      local_constraints.Aineq = sparse(0,model.size);
      local_constraints.bineq = [];
      local_constraints.lb = model.lb;      
      local_constraints.ub = model.ub;
      model.lb = -Inf(model.size,1);
      model.ub = Inf(model.size,1);
      model.local_constraints =local_constraints; 
    end
    inst_param = get_instance_run_params(model, solver, r_time, false, ...
                          epsilon, max_iter);
    inst_param.racqp_run_p = get_LPTestSet_run_params(str2num(inst(ii,3)), str2num(inst(ii,2)), ...
                     r_time, epsilon, max_iter, rnd_seed, inst(ii,4));
    s = solve_instance(inst_param);
    s.name = inst(ii,1);
    solutions = [solutions,s];
  end
  if(~quiet)
    print_solutions(solutions);
  end
end

function run_params = get_LPTestSet_run_params(n_blocks, beta, max_time, epsilon,...
                                   max_iter, rnd_seed, pL)
  
  run_params = default_run_params(n_blocks, beta,max_time,epsilon,max_iter);
  %change some parameters
  run_params.rnd_seed = rnd_seed;
  %turn on verbose
  %run_params.debug = 1;
  if(strcmpi("l_bounds",pL))
    run_params.sub_solver_type = 'gurobi';
  end
  
  gurobi_params = default_gurobi_parameters;
  %change how long gurobi will spend on each subproblem
  gurobi_params.TimeLimit = 10;
  run_params.gurobi_params = gurobi_params;
end

function inst = get_instances()

  inst=["nug30", "200", "200", "l_bounds";
    "wide15", "10", "200", "none";
    "square15", "10", "200", "none";
    "long15", "10", "200", "none";
    "i_n13", "10", "200", "none";
   "16_n14", "10", "200", "none"
];

end


