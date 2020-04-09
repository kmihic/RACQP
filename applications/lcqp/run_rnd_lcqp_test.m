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

% Randomly generated problems using procedure described in Section 4.1 
%     Mihic, K., Zhu, M. and Ye, Y., 2019. 
%     "Managing Randomization in the Multi-Block Alternating Direction 
%      Method of Multipliers for Quadratic Optimization"
%     arXiv preprint arXiv:1903.01786.


function run_rnd_lcqp_test(solver,  r_time, epsilon, max_iter, rnd_seed)

  solutions = [];
  sol_stat = [];
  inst = get_instances_lcqp();  
  dir = "../data/data_rnd/";

  for jj = 1:size(inst,1)
    sol = [];
    for ii = 1:2%10
      file_n = inst(jj)+"_RND"+ii;
      filename = file_n+".rac";
      disp("Solving: "+file_n)
      disp("LOADING THE MODEL...")
      load(dir  +file_n);
      inst_param = get_instance_run_params(model, solver, r_time, false, ...
                           epsilon, max_iter);
      inst_param.racqp_run_p = get_lcqp_run_params(str2num(inst(jj,3)), str2num(inst(jj,2)), ...
                      r_time, epsilon, max_iter, rnd_seed, model);
      s = solve_instance(inst_param);
      name = split(inst(jj),'_');
      s.name = name(2)+"_"+name(5)+"_RND"+ii;  
      sol = [sol,s];  
    end      
    solutions = [solutions,sol];  
    s = get_stats(sol);
    name = split(inst(jj),'_');
    s.name = name(2)+"_"+name(5);
    sol_stat = [sol_stat,s];
  end  
  disp(" ")
  disp("####################")
  disp(" INFO PER PROBLEM ")
  print_solutions(solutions, false, false, 'Ncol_Mrows');
  disp(" ")
  disp("####################")
  disp(" STATISTICS PER PROBLEM TYPE ")
  disp("####################")
  print_stats(sol_stat, false, false, 'Ncol_Mrows');
end


function run_params = get_lcqp_run_params(n_blocks, beta, max_time, epsilon,...
                                   max_iter, rnd_seed, model)

  run_params = default_run_params(n_blocks, beta,max_time,epsilon,max_iter);
  %change some parameters
  run_params.rnd_seed = rnd_seed;
  %turn on verbose
  %run_params.debug = 1;
  %get density
  density_Q = nnz(model.Q)/model.size^2;
  ns = size(model.Aeq,1)*size(model.Aeq,2);
  density_Aeq = nnz(model.Aeq)/ns;
  ns = size(model.Aineq,1)*size(model.Aineq,2);
  density_Aineq = nnz(model.Aineq)/ns;
  density = max([density_Q, density_Aineq, density_Aeq]);
  if(density > 0.5)
    run_params.use_sparse = false; %default value is true
  end
  
end

function inst = get_instances_lcqp()

  inst=["LCQP_N6000_SP0_95_M600_E0_5_Z2_C100", "1", "100";
  "LCQP_N6000_SP0_95_M3000_E0_5_Z2_C100", "1", "100";
  "LCQP_N9000_SP0_95_M900_E0_5_Z2_C100", "1", "150";
  "LCQP_N9000_SP0_95_M4500_E0_5_Z2_C100", "1", "150"];
end
