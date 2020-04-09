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


function run_rnd_markowitz_test(var_epsilon,r_time,  max_iter, rnd_seed)

solver = 'racqp';

  inst = get_instances();  
  solutions = [];
  sol_stat = [];
  dir = "../data/data_rnd/";

  if(var_epsilon)
    epsilon = [1e-4 1e-5 1e-6 1e-7];
    groups = 100
  else
    epsilon = 1e-5;
    groups = [50 100 150 200];    
  end
  for n_grp = groups   
    for eps = epsilon
      disp("##### NUM GROUPS: "+n_grp);
      disp("##### EPSILON: "+eps);
      sol = [];
      for jj = 1:size(inst,1)
        for ii = 1:10
          file_n = inst(jj)+"_RND"+ii;
          filename = file_n+".rac";
          disp("Solving: "+file_n)
          disp("LOADING THE MODEL...")
          load(dir+file_n);
          inst_param = get_instance_run_params(model, 'racqp', r_time, false, ...
                            eps, max_iter);
          inst_param.racqp_run_p = get_markowitz_run_params(n_grp, str2num(inst(jj,2)), ...
                       r_time, eps, max_iter, rnd_seed, model.Q);
          s = solve_instance(inst_param);
          name = split(inst(jj),'_');
          s.name = name(2)+"_RND"+ii+"_GRP_"+n_grp+"_EPS_"+eps;        
          sol = [sol,s]; 
        end %for each rnd file
      end %for each file
      solutions = [solutions,sol];  
      s = get_stats(sol);
      name = split(inst(jj),'_');
      s.name = name(2)+"_GRP_"+n_grp+"_EPS_"+eps;
      sol_stat = [sol_stat,s];
    end %for each eps     
  end  %for each n grp
  disp(" ")
  disp("####################")
  disp(" INFO PER PROBLEM ")
  print_solutions(solutions, false, false, 'Ncol_numGrp_epsilon');     
  disp(" ")
  disp("####################")
  disp(" STATISTICS PER PROBLEM TYPE ")
  disp("####################")    
  print_stats(sol_stat, false, false, 'Ncol_numGrp_epsilon');  
end

function run_params = get_markowitz_run_params(n_blocks, beta, max_time, epsilon,...
                                   max_iter, rnd_seed,Q )

  run_params = default_run_params(n_blocks, beta,max_time,epsilon,max_iter);
  %change some parameters
  run_params.rnd_seed = rnd_seed;
  %turn on verbose
  %run_params.debug = 1;
  %get density
  density = nnz(Q)/size(Q,1)^2;
  if(density > 0.5)
    run_params.use_sparse = false; %default value is true
  end

end


function inst = get_instances()
inst=["MARKOWITZ_N9000_SP0_95_E0_5_Z2_C100", "10"];

end
