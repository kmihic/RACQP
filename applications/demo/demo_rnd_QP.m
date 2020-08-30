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

% DEMO: randomly generated QP problem
%

function [solution, model] = demo_rnd_QP(problem, run_p, debug,  solver_mode)

  addpath('../../solver/racqp');
  addpath('../../solver/utils');
  addpath('../utils');

  if(nargin <= 3)
    solver_mode = 'RAC';
  end

  if(strcmpi(problem.type,'construct'))
    %construct a problem 
    model = get_model_rnd_QP(problem.construct,true);% print out the model struct
    disp(" ")
    disp("Model struct:")
    disp(model)
  elseif(strcmpi(problem.type,'model'))
    model = problem.model;
  else
    error("Problem type not correct: 'construct' or 'model' accepted only")
  end

  % clean model constraints and remove embedded slacks 
  % (handeled by RACQP as a special block)
  % No need ot run it on the generated data....
  %model = clean_model_constraints(model, true);  
  % scale model
  model = scale_model(model, true, false);

  %load default runtime parameters
  run_params = default_run_params(run_p.n_blocks, run_p.beta, run_p.max_rtime, ...
        run_p.epsilon, run_p.max_iter);
  %change some parameters
  %turn on/off verbose
  run_params.debug = debug;


  if(strcmpi('RP_ADMM', solver_mode) || strcmpi('CYCLIC_ADMM', solver_mode))
    rng(run_params.rnd_seed);
    groups={};
    g_size = floor(model.size/run_params.n_blocks);
    ix_x = randperm(model.size);
    for ii = 0:(run_params.n_blocks-1)
      st = (ii * g_size);
      en = st + g_size;
      ix = ix_x((st+1):en);
      groups((ii+1),:) = {ix, []};
    end  
    model.groups = groups;    
    if(strcmpi('RP_ADMM', solver_mode))
      model.group_mode = 'RP';
      msg = " RP-ADMM ";
    else
      model.group_mode = 'CADMM';
      msg = " CYCLIC-ADMM ";
    end
  elseif(strcmpi('RAC', solver_mode))
     msg = " RAC-ADMM ";
  else
    error("Solver mode incorrect. Accepted: RAC, RP_ADMM, CYCLIC_ADMM")
  end

  if(run_params.n_blocks == 1)
    msg = "Solver mode ignored. Running RACQP single-block";
  else
    msg = "Running RACQP multi-block:" + msg +" mode.";
  end

  %run the solver
  fprintf("\n##############\n\n%s\n\n",msg)
  solution = RACQP(model,run_params,true);

end
