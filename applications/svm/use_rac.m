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


function [r_val] = use_rac(model_libsvm,rp, no_test)
  global X_train;
  global X_train_sq;
  global Q_cached;
  global is_cached;
  global N_CACHED;
  global N_x;

  %load default runtime parameters
  N = length(model_libsvm.train_label);
  n_blocks = ceil(N/rp.max_block_size);
  beta = round(rp.beta*n_blocks,1);
  run_params = default_run_params(n_blocks,beta,rp.r_time,rp.epsilon,rp.max_iter);
  %%% change some parameters
  % Gaussian Kernel makes the model dense
  % RACQP by default does sparse calculations. 
  % Forcing the solver to do dense numeric 
  run_params.use_sparse = false;
  run_params.use_dual_res = false;
 % run_params.do_rel_tol = false;
  run_params.epsilon_dual = 1; %anything below 1 is fine!
  % calculate final dual residual (cannot be done for cached Q, not implemented Q(:,:))
  run_params.calc_dual_res = false; 
  %turn on verbose
  if(rp.verbose)
    run_params.debug = 1;
  end
  % get the model  
  disp("Preparing Data");
  [model K_train] = get_rac_model(model_libsvm, run_params, rp);

  %run the solver
  if(n_blocks == 1)
    msg = "Running RACQP single-block";
  else
    msg = "Running RACQP multi-block";
  end
  disp("Num blocks: "+n_blocks);
  disp("Beta: "+rp.beta);
  disp(msg)
  sol = RACQP(model,run_params,true); 
  % Done with optimization
  % Clear most of global vars, we are not using them any more
  clear global X_train;
  clear global X_train_sq;
  clear global Q_cached;
  clear global is_cached;

  if(no_test)
    accuracy = -1;
    disp("No testing done");
  else
    disp("Runtime: " + sol.runtime +" (RAC)")
    disp("Testing the accuracy of the training...")
    % get bias    
    alpha = round(sol.sol_x,3);
    y_alpha = alpha .* model_libsvm.train_label;    
    
    sv_ix = find(alpha>0 & alpha<model_libsvm.C); %this should be the correct way
    b = get_bias(model_libsvm.train_data(sv_ix,:),y_alpha(sv_ix),...
           model_libsvm.train_label(sv_ix));  
%    disp("Neg bias")
%    label = classify(model_libsvm.test_data, model_libsvm.train_data, y_alpha, -b);
%    correct = find(~(label - model_libsvm.test_label));
%    nc = length(correct);
%    na = length(model_libsvm.test_label);
%    accuracy = 100*nc/na;
%    msg = sprintf('Accuracy test data = %.4f%% (%d/%d) (RAC, beta=%f)',accuracy,nc, na,rp.beta);
%    disp(msg)
%    disp("Pos bias")
    label = classify(model_libsvm.test_data, model_libsvm.train_data, y_alpha, b);
    correct = find(~(label - model_libsvm.test_label));
    nc = length(correct);
    na = length(model_libsvm.test_label);
    accuracy = 100*nc/na;
    msg = sprintf('Accuracy test data = %.4f%% (%d/%d) (RAC, beta=%f)',accuracy,nc, na,rp.beta);
    disp(msg)
  end
  if(rp.static_Q)
    disp("Runtime (data preparation): " + model.data_prepare_time +" (RAC)")
    disp("Runtime (optimization): " + sol.runtime +" (RAC)")
    disp("Runtime: " + (model.data_prepare_time+sol.runtime) +" (RAC)")    
    r_time.runtime = model.data_prepare_time+sol.runtime;
    r_time.data_prepare_time = model.data_prepare_time;
    r_time.solve_time = sol.runtime;
  else
    disp("Runtime: " + sol.runtime +" (RAC)")
    r_time.runtime = sol.runtime;
    model.data_prepare_time = 0;  
    r_time.solve_time = sol.runtime;  
  end
  disp(" ")
  r_val.accuracy = accuracy;
  r_val.r_time = r_time;
%disp(" ")

end


function [model K_train] = get_rac_model(model_libsvm, run_params, rp)
global Y_train;
global X_train; %need
global YX_train; %need
global X_train_sq;
global g_gamma; 
global N_x;
  
 N_x = length(model_libsvm.train_label);
  Y_train = model_libsvm.train_label;
  train_data = model_libsvm.train_data*sqrt(g_gamma);
  X_train_sq = sum(train_data.^2,2);
  X_train = train_data;

  % Precomputed Gaussian Kernel, train
  K_train = [];
  if(rp.static_Q)
    tstart = tic;
    K_train = mygauss(model_libsvm.train_data);
    Q = bsxfun(@times, Y_train, K_train);
    model.Q = bsxfun(@times, Q, Y_train');

    t_q = toc(tstart);  
    model.data_prepare_time = t_q;
    disp("Kernel density: "+(nnz(Q)/size(Q,1)^2));
  elseif(rp.cached_Q)
    model.Q = @get_Q_cached;
  else
    model.Q = @get_Q;
  end

  N = length(model_libsvm.train_label);
  model.ub = ones(N,1)*model_libsvm.C;
  Aeq = model_libsvm.train_label';
  c = -ones(N,1);  
  model.lb = zeros(N,1);
  

  beq = 0;
  model.size = N;
  % experimenting with sparse/dense numeric
  % SVM with gaussian is a dense problem! 
  model.c = c;
  model.Aeq = Aeq;
  model.beq = 0;
  model.Aineq = zeros(0,N);
  model.bineq = [];
  % do not use any other x0 but zero. Otherwise Qx get 
  % initialized and that requires whole kernel to be processed
  model.x0 = zeros(N,1);%ones(N,1)*m.C;
  model.integers = [];
  model.binary = [];
  model.const = 0;
  model.norm_coeff = 1;

  if(strcmpi('RP', rp.racqp_mode) || strcmpi('CADMM', rp.racqp_mode))
    rng(run_params.rnd_seed);
     
    model.groups =  construct_rnd_groups(model.size, run_params.n_blocks);
    if(strcmpi('RP', rp.racqp_mode))
      model.group_mode = 'RP';
    else
      model.group_mode = 'CADMM';
    end
  end
  
end



