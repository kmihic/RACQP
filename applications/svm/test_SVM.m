function test_SVM(filename, filename_test,dir, p_train, p_test, beta, ...
        C,S, Q_model, do_rac, do_libsvm, do_matlab, solver_mode, ...
       q_proc, no_test, max_block_size, min_iter, max_iter)

  % clear global vars, just in case we did not exit properly
  clear global;
  global g_sigma;
  global g_gamma;
  global static_Q;
  global cached_Q;
  global P_CACHE;
  global USE_SWAP;

  VJENCE = false;
  USE_SWAP = false;

  if(strcmpi(Q_model,'STATIC'))
    static_Q = true;
    cached_Q = false;
  elseif(strcmpi(Q_model,'CACHED'))   
    static_Q = false;
    cached_Q = true;
    % hard coding to use 70% free memory
    P_CACHE = 0.7;
  elseif(strcmpi(Q_model,'DYNAMIC'))   
    static_Q = false;
    cached_Q = false;
  else
    error("Q model can be: STATIC, CACHED or DYNAMIC")
  end


  addpath('../RACQP/racqp');
  addpath('../LIBSVM/libsvm-master/matlab');

  if(nargin < 13)
   solver_mode = 'RAC';
  end
  if(nargin < 15)
   no_test = false;
  end
  if(nargin < 16)
   max_block_size = 100;
  end
  if(nargin < 17)
   min_iter = 2;
   max_iter = 10;
  end
  if(nargin < 14)
    q_proc = false;
  end
 

  

  rng(1); 
  %set matlab threads to 1
  maxNumCompThreads(1);

  disp(" ")
  disp("Solving "+filename);

  epsilon = 1e-1;
  acc_RAC = [];
  time_RAC = [];
  sigma_RAC = [];
  c_RAC = [];
  acc_LIB = [];
  time_LIB = [];
  acc_MAT = [];
  time_MAT = [];
  beta_RAC = [];
  S_run = [];
  C_run = [];
  ii = 1;
  for sigma = S; 
    g_sigma = sigma;
    g_gamma = 1/(2*sigma^2);
    % prapare data
    file = [dir,filename];
    if(~strcmpi(filename_test,''))
      fileT = [dir,filename_test];
    else
      fileT="";
    end
    m = prepare_data(file, p_train, fileT, p_test); 
    m.sigma = sigma;
    % do different C
    [sol_RAC, sol_LIB, sol_MAT] = solve_SVM(m, C, epsilon, beta, do_rac, ...
       do_libsvm, do_matlab, solver_mode, no_test,max_block_size, min_iter, max_iter );    
    C_run = [C_run C];
    S_run = [S_run ones(1,length(C))*sigma];
    if(do_rac)
      acc_RAC = [acc_RAC sol_RAC.accuracy];
      time_RAC = [time_RAC sol_RAC.r_time];
      beta_RAC = [beta_RAC sol_RAC.beta];
      sigma_RAC = [sigma_RAC sol_RAC.sigma];
      c_RAC = [c_RAC sol_RAC.C];
    end
    if(do_libsvm)
      acc_LIB = [acc_LIB sol_LIB.accuracy];
      time_LIB = [time_LIB sol_LIB.r_time];
    end
    if(do_matlab)
      acc_MAT = [acc_MAT sol_MAT.accuracy];
      time_MAT = [time_MAT sol_MAT.r_time];
    end
  end

  %T = table(name,obj_val,run_time);
  %T.Properties.VariableNames={'Data_file','obj_val', 'run_time'};
  %disp(T);
  if(do_rac)
    disp(solver_mode)
    T_runp = table(beta_RAC',c_RAC', sigma_RAC', acc_RAC');
    T_runp.Properties.VariableNames={'Beta', 'C', 'Sigma','Accuracy'};
    T_runt=struct2table(time_RAC);
    T = [T_runp T_runt];
    disp(T)
  end
  if(do_libsvm)
    disp('LIBSVM')
    T = table(C_run', S_run',acc_LIB', time_LIB');
    T.Properties.VariableNames={'C', 'Sigma','Accuracy','libsvm_time'};
    disp(T)    
  end

  if(do_matlab)
    disp('MATLAB')
    T = table(C_run', S_run',acc_MAT',time_MAT');
    T.Properties.VariableNames={'C', 'Sigma','Accuracy','matlab_time'};
    disp(T)    
  end
  if(q_proc)
    quit
  end
  clear all;
end


function [sol_RAC, sol_LIBSVM, sol_MAT] = solve_SVM(m, C, epsilon, beta, ...
       do_rac, do_libsvm, do_matlab, solver_mode, no_test, max_block_size,...
       min_iter, max_iter)

  %run RAC    
  global g_sigma;
    N = length(m.train_label);
    block_size = max_block_size;
    n_blocks = ceil(N/block_size);

    
  r.n_blocks = n_blocks;
  r.epsilon = epsilon;
  r.max_iter = max_iter;
  r.min_iter = min_iter;
  r.max_time = 360000;
  r.verbose = true;
  
  disp("Train set size: "+N);
  disp("Test set size: "+length(m.test_label));
  disp("Num features: "+size(m.train_data,2));
  disp(" ")
  ii = 1;
  sol_RAC = [];
  sol_LIBSVM = [];
  sol_MAT=[];
  for c = C   
    msg = sprintf("\n------ sigma = %f, c = %f --------\n",m.sigma,c);
    disp(msg);
    m.C = c;
    % run rac
    if(do_rac)
     for b = beta
       if(beta<0)
        r.beta = -beta;
       else
        r.beta = round(b*n_blocks,1);
       end
       s = use_rac(m,r, solver_mode, no_test);
       s.beta = b;
       s.C = c;
       s.sigma = g_sigma;
       sol_RAC=[sol_RAC s];
     end
    end
    if(do_libsvm)
      %run libsvm    
      sol_LIBSVM = [sol_LIBSVM use_libsvm(m, no_test)];
    end
    if(do_matlab)
      %run matlab
      sol_MAT = [sol_MAT use_matlab(m, no_test)];
    end
  end
  
end
