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


function run_LPTestSet_test()
%addpath('/opt/gurobi/gurobi752/linux64/matlab/');
addpath('../racqp');
addpath('../utils');

  inst = get_instances();
  solutions = [];
  for ii = 1:length(inst)
    disp('Solving: '+inst(ii,1));

    s = verify_LPTestSet(inst(ii,1),str2num(inst(ii,2)),str2num(inst(ii,3)),inst(ii,4));
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
  
  
%quit
end

function inst = get_instances()
inst=["../data/data_mps/wide15.mps", "10", "200", "none";
"../data/data_mps/square15.mps", "10", "200", "none";
"../data/data_mps/long15.mps", "10", "200", "none";
"../data/data_mps/nug30.mps", "200", "200", "l_bounds"];
end


