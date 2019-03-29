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


function run_min_variance_test()
addpath('../racqp');
addpath('../utils');


  inst = get_instances();
  solutions = [];
  for ii = 1:length(inst)
    disp('Solving: '+inst(ii,1));
    s = verify_Markowitz_portfolio_selection_data(inst(ii,1),inst(ii,2),str2num(inst(ii,3)));
    solutions = [solutions,s];
  end
  disp(" ")
  disp("#####################")
  disp('SUMMARY')
  for ii = 1: length(solutions)
    name = inst(ii,1);
    obj_val = solutions(ii).sol_obj_val;
    run_time = solutions(ii).rac_time;
    out = name + " " + obj_val + " " + run_time;
    disp(out)
  end
end


function inst = get_instances()
%regular files are to large for gitHub
%inst=["../data/data_markowitz/regular_daily.mat","false","60";  
%"../data/data_markowitz/regular_monthly.mat","false","60";
%"../data/data_markowitz/regular_quarterly.mat","false","60";
%"../data/data_markowitz/lowrank_daily.mat","true","30";  
%"../data/data_markowitz/lowrank_monthly.mat","true","30";  
%"../data/data_markowitz/lowrank_quarterly.mat","true","30"];
inst=["../data/data_markowitz/lowrank_daily.mat","true","30";  
"../data/data_markowitz/lowrank_monthly.mat","true","30";  
"../data/data_markowitz/lowrank_quarterly.mat","true","30"];
end
