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

function run_svn_test(solver)

  addpath('../../solver/racqp');
  addpath('../../solver/utils');
  addpath('../../libsvm/libsvm-master/matlab');

  % clear global vars, just in case we did not exit properly
  global g_gamma;

  rng(1); 
  %set matlab threads to 1
  maxNumCompThreads(1);

  no_test = false;
  dir = "../data/data_libsvm/";
  Q_model = 'cached';
  p_train = 1;
  p_test = 1;
  racqp_mode = 'RAC';
  beta = 0.1;

  if(~(strcmpi(solver,'racqp') ...
     || strcmpi(solver,'libsvm') ...
     || strcmpi(solver,'matlab')))
    error("Solver not recognised. Allowed: 'racqp', libsvm' and 'matlab'");
  end
  rac_rp = prepare_rac_params(Q_model, racqp_mode,beta);
  inst = get_instances_svm();
  solutions = [];
  for ii = 1:size(inst,1)
    disp(" ")
    disp("Solving "+inst(ii,1));
    disp("Loading the model");
    file_n = dir+inst(ii,1);
    load(file_n);
    sigma = str2num(inst(ii,2));
    c = str2num(inst(ii,3));
    rac_rp.max_block_size = str2num(inst(ii,4));
    do_matlab = str2num(inst(ii,5));
    disp("Train set size: "+length(model_libsvm.train_label));
    disp("Test set size: "+length(model_libsvm.test_label));
    disp("Num features: "+size(model_libsvm.train_data,2));
    disp(" ");
    if(~do_matlab && strcmpi(solver,'matlab'))
      disp("Instance too large for Matlab or more than 10h needed to solve it")
      disp("Skipping the instance");
    else 
      g_gamma = 1/(2*sigma^2);    
      model_libsvm.sigma = sigma;
      model_libsvm.C = c;
      sol = solve_instance(model_libsvm, solver, rac_rp, no_test,beta);
      sol.name = inst(ii,1);
      solutions = [solutions, sol];     
    end   
  end
  print_solutions_svm(solutions);
  clear global;
quit
end

function rac_rp = prepare_rac_params(Q_model, racqp_mode,beta)

% racqp params
  rac_rp = struct();
  rac_rp.verbose = false;
  rac_rp.r_time = Inf;
  rac_rp.max_iter = 5;
  rac_rp.epsilon = 1e-1;
  rac_rp.racqp_mode = racqp_mode;
  
  if(strcmpi(Q_model,'static'))
    rac_rp.static_Q = true;
    rac_rp.cached_Q = false;
  elseif(strcmpi(Q_model,'cached'))   
    rac_rp.static_Q = false;
    rac_rp.cached_Q = true;
  elseif(strcmpi(Q_model,'dynamic'))   
    rac_rp.static_Q = false;
    rac_rp.cached_Q = false;
  else
    error("Q model can be: static, cached or dynamic.")
  end 


end

function inst = get_instances_svm()
inst=[
"a8a","0.1","1","100", "1"; 
"w7a","1","0.1","100", "1"; 
"rcv1_binary","10","0.1","100", "0"; 
"news20_binary","0.1","1","100", "0"; 
"a9a","0.1","1","500", "1"; 
"w8a","0.1","1","500", "1"; 
"ijcnn1","0.1","0.1","500", "1"; 
"cod_rna","0.1","0.1","500", "1";  
"real_sim","0.1","0.1","500", "0"; 
"skin_nonskin","10","0.1","500", "0";
"webspam_uni","10","10","100","0"; 
"covtype_binary","0.1","10","1000", "0";
];
end

