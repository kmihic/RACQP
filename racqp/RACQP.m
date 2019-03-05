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

function rac_out = RACQP(model, run_p)
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

  disp("This version of RACQP implements")
  disp("   min  x'H'x + c'x")
  disp("   s.t. Aeq x = beq")
  disp("        Aineq x <= bineq")
  disp("        lb <= x <= ub")
  disp("        x integer, continuous")
  disp("")
  %set random seed for this run
  rng(run_p.rnd_seed);
  %set matlab threads to 1
  maxNumCompThreads(1);

  %check for probem type
  use_mip = false;
  if((length(model.integers) > 0 || length(model.binary) > 0))
    use_mip = true;
  end
  
  %remove embedded slacks to avoid problem with cholesky
  if(~use_mip && strcmpi('cholesky',run_p.sub_solver_type))
    model = remove_embedded_slacks(model);
  end
  %call the solver
  if(use_mip)
    run_p.run_sub.mip=use_mip;
    rac_out = rac_mip(model, run_p);
  elseif(run_p.n_blocks == 1)
    run_p.mip=false;
    rac_out = rac_single_block(model, run_p);
  else
    run_p.mip=false;
    rac_out = rac_multi_block(model, run_p);
  end
  
  %adjust the objective value, if objective scaled
  if(isfield(model,'norm_coeff'))
    rac_out.sol_obj_val = rac_out.sol_obj_val*model.norm_coeff;
  end
end

