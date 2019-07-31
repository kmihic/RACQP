
function [r_val] = use_rac(m,r, solver_mode, no_test)
  global static_Q;
  global X_train;
  global X_train_sq;
  global Q_cached;
  global is_cached;
  global N_CACHED;
  global N_x;

  %load default runtime parameters
  run_params = default_run_params;
  %%% change some parameters
  % Gaussian Kernel makes the model dense
  % RACQP by default does sparse calculations. 
  % Forcing the solver to do dense numeric 
  M_SPARSE = false;
  run_params.use_sparse = M_SPARSE;

  run_params.n_blocks = r.n_blocks;
  run_params.beta = r.beta;
  run_params.epsilon = r.epsilon;
  %default min number of iterations is >2
  %if less, we may get lucky and hit primal residual less then the tolerance
  %but we may be far far away from the optimal pt. We should have at least 
  %a couple of runs
  run_params.max_iter = r.max_iter;
  run_params.min_iter = r.min_iter;
  % do not claculate dual residual. Anyway it is used for output info only
  % and takes away time
  run_params.no_calc_dual_res = true;
  run_params.max_rtime = r.max_time;
  run_params.clean_embedd = false;
  %turn off verbose
  if(r.verbose)
    run_params.debug = 2;
  else
    run_params.debug = 0;
  end

  if(r.n_blocks == 1)
    msg = "Running RACQP single-block";
  else
    msg = "Running RACQP multi-block:";
  end
  disp("Num blocks: "+r.n_blocks);
  disp("Beta: "+r.beta);

  % get the model
  model = get_rac_model(m, M_SPARSE, run_params, solver_mode);

%model.group_mode = 'RP';
%model.groups=[];
  %run the solver
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
    disp("Runtime: " + sol.rac_time +" (RAC)")
    % get bias    
    alpha = sol.sol_x.* m.train_label;
    alpha(alpha<0) = 0;
    y_alpha = alpha .* m.train_label;
    
    sv_ix = find(alpha>0); 
    sv_c = find(alpha>0 & alpha<m.C);
    %b = get_bias(m.train_data(sv_c,:),y_alpha(sv_c), m.train_label);
    b = get_bias(m.train_data(sv_ix,:),y_alpha(sv_ix), m.train_label);
    
    label = classify(m.test_data, m.train_data(sv_ix,:), y_alpha(sv_ix), b);
    %label = classify(m.test_data, m.train_data, y_alpha, b);

    correct = find(~(label - m.test_label));
    nc = length(correct);
    na = length(m.test_label);
    accuracy = 100*nc/na;
  
    disp("RAC test data")
    msg = sprintf('Accuracy test data = %.4f%% (%d/%d) (RAC, beta=%f)',accuracy,nc, na,r.beta);
    disp(msg)
  end
  if(static_Q)
    disp("Runtime (data preparation): " + model.data_prepare_time +" (RAC)")
    disp("Runtime (optimization): " + sol.rac_time +" (RAC)")
    disp("Runtime: " + (model.data_prepare_time+sol.rac_time) +" (RAC)")
  else
    disp("Runtime: " + sol.rac_time +" (RAC)")
    model.data_prepare_time = 0;
  end
  disp(" ")
% disp("RP factorization: "+sol.rp_time_factor)
% disp("Getting Q(:,x_ix): "+sol.r_time_Q)
% disp("Calc Q(:,x_ix)*x: "+sol.r_time_Qx)
% disp("Sub-problem solving: "+sol.solver_time)


  if(static_Q)
    r_time.Q_static = model.data_prepare_time;
  end
%   r_time.rp_time_factor = sol.rp_time_factor;
%   r_time.r_time_Q = sol.r_time_Q;
%   r_time.r_time_Qx = sol.r_time_Qx;
%   r_time.sub_solver_time = sol.solver_time;
   r_time.rac_time = sol.rac_time;  
  if(static_Q)
    r_time.total = model.data_prepare_time+sol.rac_time;
  end

r_val.accuracy = accuracy;
r_val.r_time = r_time;
disp(" ")

end


function model = get_rac_model(m, m_sparse, run_params, solver_mode)
%global OLD;
%global QY;
global Y_train;
global X_train; %need
global YX_train; %need
global X_train_sq;
global static_Q; %need
global cached_Q;
global g_gamma; 
global N_x;
%global YY;
  
 N_x = length(m.train_label);
  Y_train = m.train_label;
  train_data = m.train_data*sqrt(g_gamma);
  X_train_sq = sum(train_data.^2,2);
  if(~m_sparse)
    X_train_sq = full(X_train_sq);
  end
  X_train = train_data;

  % Precomputed Gaussian Kernel, train
  if(static_Q)
    tstart = tic;
    K_train = mygauss2(m.train_data);
    Q = bsxfun(@times, Y_train, K_train);
    model.Q = bsxfun(@times, Q, Y_train');
    t_q = toc(tstart);  
    model.data_prepare_time = t_q;
    disp("Runtime (data preparation): " + t_q +" (RAC)"); 
    disp("Kernel density: "+(nnz(Q)/size(Q,1)^2));
  elseif(cached_Q)
    model.Q = @get_Q_cached;
  else
    model.Q = @get_Q;
  end

  N = length(m.train_label);
  model.ub = ones(N,1)*m.C;
    c = -m.train_label;
    Aeq = ones(1,N);  
    model.lb = -model.ub;
%   end
  

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

  if(strcmpi('RP_ADMM', solver_mode) || strcmpi('CYCLIC_ADMM', solver_mode))
    rng(run_params.rnd_seed);
     
    model.groups = make_groups(run_params.n_blocks,model.size);    
    if(strcmpi('RP_ADMM', solver_mode))
      model.group_mode = 'RP';
    else
      model.group_mode = 'CADMM';
    end
  end
  
end



