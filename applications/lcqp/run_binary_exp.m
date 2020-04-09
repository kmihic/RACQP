%
% RACQP -  Randomly Assembled Cyclic ADMM Quadratic Programming Solver
% Copyright (C) 2019 
%     Kresimir Mihic <kmihic@alumni.stanford.edu>
%     Mingxi Zhu <mingxiz@stanford.edu>
%     Yinyu Ye <yyye@stanford.edu>
%
% This file is part of RACQP 
%

% Verify binary results 



function run_binary_exp(section,solver, r_time)
 
  addpath('/opt/mosek/8/toolbox/r2014a');
  addpath('/opt/gurobi/gurobi811/linux64/matlab/');
  addpath('/opt/osqp-matlab');
  addpath('../../solver/racqp');
  addpath('../../solver/utils');
  addpath('../utils');
  warning off;

  if(nargin <= 2)
    r_time = -1;
  end
  if(r_time > 0 || (r_time == -1 && strcmp(solver,'gurobi')) || r_time == 10 ...
      || r_time == 60 || r_time == 300 || r_time == 600 ...
      || r_time == 1800 || r_time == 3600)        
    if(strcmp(solver,'racqp') || strcmp(solver,'gurobi') )
      solve(section, solver, r_time);
    else
      disp("Solver not recognized");
      disp("Accepted: 'racqp', 'gurobi'");
    end
   else
    disp("Incorrect call")
    disp("Allowed values for r_time: 60, 300, 600 and 1800 (s)");
   end 
 

  quit
end

function solve(section, solver, r_time)

  binary=true;
  switch section
    case '4.2.2_regular'
      msg = "4.2.2 Markowitz Portfolio Selection, regular";
      disp("Solving section "+msg);
      disp("Solver: "+solver);
      disp(" ");        
      run_portfolio_test(solver, r_time, false, binary);
    case '4.2.2_lowrank'
      msg = "4.2.2 Markowitz Portfolio Selection, low rank";
      disp("Solving section "+msg);
      disp("Solver: "+solver);
      disp(" ");        
      run_portfolio_test(solver, r_time, true, binary); 
    case '4.2.3_full'  
      msg = "4.2.3 QAPLIB binary, full benchmark";
      disp("Solving section "+msg);
      disp("Solver: "+solver);
      disp(" ");        
      run_qaplib_test(solver,r_time,binary,true);
    case '4.2.3_large'  
      msg = "4.2.3 QAPLIB binary, large instances";
      disp("Solving section "+msg);
      disp("Solver: "+solver);
      disp(" ");        
      run_qaplib_test(solver,r_time,binary,false);
    case '4.2.4'
      msg = "4.2.4 Maximum Cut Problem";
      disp("Solving section "+msg);
      disp("Solver: "+solver);
      disp(" ");
      run_Gcut_test(solver, r_time, false);
    case '4.2.5'
      msg = "4.2.5 Maximum Bisection Problem";
      disp("Solving section "+msg);
      disp("Solver: "+solver);
      disp(" ");
      run_Gcut_test(solver, r_time, true);
    otherwise
      disp("Section not recognized");
      quit
  end
end
