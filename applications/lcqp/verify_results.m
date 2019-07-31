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

% Verify results using RACQP

function verify_results(section,solver, r_time)
 
  % Note: For some continious instances and all binary problems RACQP requires gurobi
  %addpath('/opt/gurobi/gurobi752/linux64/matlab/');
  %addpath('/opt/mosek/8/toolbox/r2014a');
  if(nargin == 2)
    r_time = -1;
  end
  if(strcmp(solver,'racqp') || strcmp(solver,'gurobi') || strcmp(solver,'mosek'))
    solve(section, solver, r_time);
  else
    disp("Solver not recognized");
    disp("Accepted: 'racqp', 'gurobi', 'mosek'");
  end

  %quit
end

function solve(test, solver, r_time)
  
  switch test
    case '4.1.1'
      msg = "4.1.1 Regularized Markowitz Mean-Variance Model";
      disp("Solving section "+msg+" using RACQP");
      disp(" ");
      run_min_variance_test(solver)
    case '4.1.2'
      disp("Must define subsection")
      disp("Allowed: 4.1.2_markow, 4.1.2_lcqp")
    case '4.1.2_markow'
      msg = "4.1.2 Randomly Generated Linearly Constrained Quadratic Problems (LCQP)";      
      disp("Solving section "+msg+" using RACQP");
      disp("Markowitz-like Problem Instances");
      disp(" ");
      run_rnd_markowitz_like_test(solver)
    case '4.1.2_lcqp'
      msg = "4.1.2 Randomly Generated Linearly Constrained Quadratic Problems (LCQP)";      
      disp("Solving section "+msg+" using RACQP");
      disp("General LCQP")
      disp(" ");
      run_rnd_LCQP_test(solver)
    case '4.1.3'
      msg = "4.1.3 Relaxed QAP";
      disp("Solving section "+msg+" using RACQP");
      disp(" ");
      run_qap_relaxed_test(solver)
    case '4.1.4'
      msg = "4.1.4 Maros and Meszaros Convex QP";
      disp("Solving section "+msg+" using RACQP");
      disp(" ");
      run_cuter_test(solver)
    case '4.1.5'
      msg = "4.1.5 Convex QP based on the Mittelmann LP test set"
      disp("Solving section "+msg+" using RACQP");
      disp(" ");
      run_LPTestSet_test(solver)
    case '4.2.1'
      % 4.2.1 Randomly Generated, Linearly Constrained Binary Quadratic Problems
      disp("Data files huge. Java based generator used to make the files. Avaiable upon request.")
    case '4.2.2'
      disp("Must define subsection")
      disp("Allowed: 4.2.2_regular, 4.2.2_lowrank")
    case '4.2.2_regular'
      msg = "4.2.2 Markowitz Portfolio Selection";
      disp("Solving section "+msg+" using RACQP");
      disp(" ");
      if((strcmp(solver,'racqp') || strcmp(solver,'gurobi')) ...
         && (r_time == 60 || r_time == 300 || r_time == 600 || r_time == 1800))        
        run_min_portfolio_test(solver, r_time)
      else
        disp("Incorrect call")
        disp("verify_results('4.2.2',solver,r_time)")
        disp("Allowed values for r_time: 60, 300, 600 and 1800 (s)");
        disp("Allowed values for solver: 'racqp', 'gurobi'");
      end      
    case '4.2.2_lowrank'
      msg = "4.2.2 Markowitz Portfolio Selection";
      disp("Solving section "+msg+" using RACQP");
      disp(" ");
      run_min_portfolio_lr_test(solver)
    case "4.2.3"
      msg = "4.2.3 QAPLIB";
      disp("Solving section "+msg+" using RACQP");
      disp(" ");
      if(strcmp(solver,'racqp') || strcmp(solver,'gurobi'))
        run_qap_binary_test(solver)
      else
        disp("Incorrect call")
        disp("verify_results('4.2.3',solver)")
        disp("Allowed values for solver: 'racqp', 'gurobi'");
      end
    case "4.2.4"
      msg = "4.2.4 Maximum Cut Problem";
      disp("Solving section "+msg+" using RACQP");
      disp(" ");
      if((strcmp(solver,'racqp') || strcmp(solver,'gurobi')) ...
         && (r_time == 600 || r_time == 1800 || r_time == 3600))
        run_maxcut_test(r_time,solver);
      else
        disp("Incorrect call")
        disp("verify_results('4.2.4',solver,r_time)")
        disp("Allowed values for r_time: 600 1800 and 3600 (s)");
        disp("Allowed values for solver: 'racqp', 'gurobi'");
      end
    case "4.2.5"
      msg = "4.2.5 Maximum Bisection Problem";
      disp("Solving section "+msg+" using RACQP");
      disp(" ");
      if((strcmp(solver,'racqp') || strcmp(solver,'gurobi')) ...
         && (r_time == 600 || r_time == 1800 || r_time == 3600))
        run_bisection_test(r_time,solver);
      else
        disp("Incorrect call")
        disp("verify_results('4.2.5',solver,r_time)")
        disp("Allowed values for r_time: 600 1800 and 3600 (s)");
        disp("Allowed values for solver: 'racqp', 'gurobi'");
      end
    otherwise
      error("Section not recognized");
  end
end



