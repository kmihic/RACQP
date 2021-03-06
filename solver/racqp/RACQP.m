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

function rac_out = RACQP(model, run_p, quiet)
%RACQP RAC-ADMM Solver for Quadratic Problems 
%   RACQP solves QP models of the following format
%         min  x'H'x + c'x
%         s.t. Aeq x = beq
%              Aineq x <= bineq
%              lb <= x <= ub
%              x integer, continuous
%
%   rac_out = RACQP(model, run_parameters) solves QP model using runtime parameters
%
  time_start = tic;
  if(nargin <3 || ~quiet)
    disp("This version of RACQP implements")
    disp("   min  x'H'x + c'x")
    disp("   s.t. Aeq x = beq")
    disp("        Aineq x <= bineq")
    disp("        lb <= x <= ub")
    disp("        x integer, continuous")
    disp("")
  end
  %set random seed for this run
  rng(run_p.rnd_seed);
  %set matlab threads to 1
  maxNumCompThreads(1);

  %check for probem type
  use_mip = false;
  if((length(model.integers) > 0 || length(model.binary) > 0))
    use_mip = true;
  end
  
   %if RP or ADMM modes selected and no groups given in the model
   %construct them at random
   if(isfield(model,'group_mode') ...
      && (strcmpi('RP',model.group_mode) ...
          || strcmpi('CADMM',model.group_mode) ...
          || strcmpi('DADMM',model.group_mode)) ...
      && ~isfield(model,'groups'))
    if(use_mip)
      model.groups = construct_rnd_groups(model.size,run_p.run_sub.n_blocks);
    else
      model.groups = construct_rnd_groups(model.size,run_p.n_blocks);
    end
  end
  %call the solver
  if(use_mip)
    run_p.run_sub.mip=use_mip;
    run_p.run_sub.epsilon_prim = run_p.mip_epsilon;
    run_p.run_sub.epsilon_dual = 0;
    run_p.run_sub.do_rel_tol = false;
    rac_out = rac_mip(model, run_p,time_start);
  elseif(run_p.n_blocks == 1)
    run_p.mip=false;
    rac_out = rac_single_block(model, run_p, time_start);
  else
    run_p.mip=false;
    rac_out = rac_multi_block(model, run_p, time_start);
  end
  
  %adjust the objective value, if objective scaled
  if(isfield(model,'norm_coeff'))
    rac_out.sol_obj_val = rac_out.sol_obj_val*model.norm_coeff;
  end
end

