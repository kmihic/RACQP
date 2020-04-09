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
% Verify results: Max-cut and Max-bisection problems
%      Using GSET instances 
%      (http://web.stanford.edu/ yyye/yyye/Gset)
%
%


function run_Gcut_test(solver, r_time, bisection)

  max_time = r_time;
  inst = get_instances();

  solutions = [];
  data_path = "../data/data_gset/";
  for ii = 1:size(inst,1)
    filename = data_path+inst(ii,1);
    disp('Solving: '+filename);
    disp("LOADING THE MODEL...")
    if(bisection)
      model = load_GSET(filename,123,0.5);
    else
      model = load_GSET(filename);
    end
    inst_param = get_instance_run_params(model, solver, r_time, true);
    inst_param.racqp_run_p = get_Gcut_run_params(model.size, max_time,bisection);
    s = solve_instance(inst_param);
    %s = verify_Gcut(filename, max_time,solver, bisection);
    %this is a maximization problem, flipping the obj value
    s.sol_obj_val = -s.sol_obj_val;    
    s.name = inst(ii,1);
    if(bisection)
      s.obj_val = inst(ii,3);
    else
      s.obj_val = inst(ii,2);
    end
    solutions = [solutions,s];
  end
  disp(" ")
  disp("#####################")
  disp('SUMMARY')  
  print_solutions_binary(solutions)

end

function run_mip = get_Gcut_run_params(nsize, max_time, bisection)

  % run-parameters
  n_blocks = 4;
  p_trial = 0.005;
  lambda = 0.4;
  beta = 0.005;
  
  %load default runtime parameters
  run_params = default_run_params(n_blocks, beta,max_time);
  %change some parameters
  %turn on verbose
  %run_params.debug = 1;

  %change some parameters
  run_params.rnd_seed = 123;
  if(bisection)
    run_params.sub_solver_type='user_defined';
    run_params.sub_solver_f = @solve_subproblem_gset;
  else
    run_params.sub_solver_type='gurobi';
  end
  
  %load default gurobi parameters
  gurobi_params = default_gurobi_parameters;
  run_params.gurobi_params = gurobi_params;
  
  %set MIP run parameters
  if(bisection)
    run_mip.permute_type = 'swap';
  else
    run_mip.permute_type = 'permute';
  end
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
  run_mip.mip_epsilon = 0;
end
  
  
  
  

function inst = get_instances()
inst=[% inst name, max-cut opt, mab-bisect opt
"G1","11624","11624"; 
"G6","2178","2177"; 
"G11","564","564"; 
"G14","3064","3062"; 
"G18","992","992"; 
"G22","13359","13359"; 
"G27","3848","3341"; 
"G32","1410","1410"; 
"G36","7678","7678"; 
"G39","2408","2408"; 
"G43","6660","6659"; 
"G50","5880","5880"; 
"G51","3848","3847"; 
"G55","10299","10299"; 
"G56","4016","4016"; 
"G58","19276","19276"; 
"G60","14187","14187"; 
"G61","5796","5796"; 
"G63","26997","26988"; 
"G67","6940","6938"; 
"G70","9581","9581"; 
"G77","9926","9918"; 
"G81","14030","14030";
];
end

