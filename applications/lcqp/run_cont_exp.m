%
% RACQP -  Randomly Assembled Cyclic ADMM Quadratic Programming Solver
% Copyright (C) 2019 
%     Kresimir Mihic <kmihic@alumni.stanford.edu>
%     Mingxi Zhu <mingxiz@stanford.edu>
%     Yinyu Ye <yyye@stanford.edu>
%
% This file is part of RACQP 
%

% Verify cont results 


function run_cont_exp(section,solver, r_time, epsilon, max_iter, rnd_seed)
 
  addpath('/opt/mosek/8/toolbox/r2014a');
  addpath('/opt/gurobi/gurobi811/linux64/matlab/');
  addpath('/opt/osqp-matlab');
  addpath('../../solver/racqp');
  addpath('../../solver/utils');
  addpath('../utils');
  warning off;

  if(nargin <= 2)
    r_time = 10800;
  end
  if(nargin <= 3)
    %all experiments except LPTestSet and Cute use 1e-5
    if(strcmpi(section,'4.1.5') || strcmpi(section,'4.1.6'))
      epsilon = 1e-4;
    else
      epsilon = 1e-5;
    end
  end
  if(nargin <= 4)
    max_iter = 4000;
  end
  if(nargin <= 5)
    rnd_seed = 123;
  end

  if(strcmp(solver,'racqp') || strcmp(solver,'gurobi')... 
     || strcmp(solver,'osqp') || strcmp(solver,'mosek') )
    solve(section, solver, r_time, epsilon, max_iter, rnd_seed);
  else
    disp("Solver not recognized");
    disp("Accepted: 'racqp', 'gurobi', 'osqp' and 'mosek'");
  end

%  quit
end

function solve(section, solver, r_time, epsilon, max_iter, rnd_seed);

  binary=false;
  switch section 
    case 'admm_variants'    
      if(~strcmpi(solver,'racqp'))
       disp("ADMM variants comparison can be ran by RACQP only!")
       quit
      end
      msg = "ADMM variants comparison";
      disp(msg);
      disp(" ");
      run_admm_comparison_test(Inf,rnd_seed);
    case '4.1.2_regular'
      msg = "4.1.2 Markowitz min variance, regular";
      disp("Solving section "+msg);
      disp("Solver: "+solver);
      disp(" ");       
      % special case for osqp in this section
      if(strcmpi(solver,'osqp') && max_iter < 20000)
        max_iter = 20000;
      end
      run_portfolio_test(solver, r_time, false, binary, epsilon, max_iter, rnd_seed);
    case '4.1.2_lowrank'
      msg = "4.1.2 Markowitz min variance,  low rank";
      disp("Solving section "+msg);
      disp("Solver: "+solver);
      disp(" ");        
      run_portfolio_test(solver, r_time, true, binary, epsilon, max_iter, rnd_seed); 
    case '4.1.3_markow_grps' %
      if(~strcmpi(solver,'racqp'))
       disp("Section 4.1.3 can be ran by RACQP only!")
       quit
      end
      msg = "4.1.3 Randomly Generated Linearly Constrained Quadratic Problems (Markow). Performance: num groups.";
      disp("Markowitz-like Problem Instances");
      disp("Solving section "+msg);
      disp("Solver: "+solver);
      disp(" ");
      run_rnd_markowitz_test(false,r_time,  max_iter, rnd_seed);
    case '4.1.3_markow_epsilon' %
      if(~strcmpi(solver,'racqp'))
       disp("Section 4.1.3 can be ran by RACQP only!")
       quit
      end
      msg = "4.1.3 Randomly Generated Linearly Constrained Quadratic Problems (Markow). Performance: epsilon.";
      disp("Markowitz-like Problem Instances");
      disp("Solving section "+msg);
      disp("Solver: "+solver);
      disp(" ");
      run_rnd_markowitz_test(true,r_time,  max_iter, rnd_seed);
    case '4.1.3_lcqp' %
      msg = "4.1.3 Randomly Generated Linearly Constrained Quadratic Problems (LCQP)";
      disp("General LCQP")
      disp("Solving section "+msg);
      disp("Solver: "+solver);
      disp(" ");
      run_rnd_lcqp_test(solver,  r_time, epsilon, max_iter, rnd_seed);
    case '4.1.4' %
      msg = "4.1.4 Relaxed QAP";
      disp("Solving section "+msg);
      disp("Solver: "+solver);
      disp(" ");
      run_qaplib_test(solver,r_time,binary, false, epsilon, max_iter, rnd_seed);
    case '4.1.5' %
      msg = "4.1.5 Maros and Meszaros Convex QP";
      disp("Solving section "+msg);
      disp("Solver: "+solver);
      disp(" ");
      run_cuter_test(solver, r_time, epsilon, max_iter, rnd_seed);
    case '4.1.6' %
      msg = "4.1.6 Convex QP based on the Mittelmann LP test set";
      disp("Solving section "+msg);
      disp("Solver: "+solver);
      disp(" ");
      run_LPTestSet_test(solver, r_time, epsilon, max_iter, rnd_seed);
    case '4.1.7' %
      if(~strcmpi(solver,'racqp'))
       disp("Section 4.1.7 can be ran by RACQP only!")
       quit
      end
      msg = "4.1.7 RAC rnd seed";
      disp("Solving section "+msg);
      disp("Solver: "+solver);
      disp(" ");
      run_rnd_seed_test(solver, r_time, epsilon, max_iter);
    otherwise
      disp("Section not recognized");
      quit
  end
end
