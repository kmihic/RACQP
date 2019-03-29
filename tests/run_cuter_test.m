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


function run_cuter_test()
%addpath('/opt/gurobi/gurobi752/linux64/matlab/');
addpath('../racqp');
addpath('../utils');

  inst = get_instances();
  solutions = [];
  for ii = 1:length(inst)
    disp('Solving: '+inst(ii,1));
    s = verify_CUTEr(inst(ii,1),inst(ii,2),inst(ii,3),inst(ii,4));
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
inst=["../data/data_cuter/AUG2D","m","none","false";
"../data/data_cuter/AUG2DC","m","none","false";
"../data/data_cuter/AUG2DCQP","m","l_bounds","false";
"../data/data_cuter/AUG2DQP","m","l_bounds","false";
"../data/data_cuter/AUG3D","m","none","false";
"../data/data_cuter/AUG3DC","m","none","false";
"../data/data_cuter/AUG3DCQP","m","l_bounds","false";
"../data/data_cuter/AUG3DQP","m","l_bounds","false";
"../data/data_cuter/CONT-050","m","l_bounds","false";
"../data/data_cuter/CONT-100","m","l_bounds","false";
"../data/data_cuter/CONT-101","m","l_bounds","false";
"../data/data_cuter/CONT-200","m","l_bounds","false";
"../data/data_cuter/CONT-201","m","l_bounds","false";
"../data/data_cuter/CONT-300","m","l_bounds","false";
"../data/data_cuter/CVXQP1_L","m3","l_bounds","false";
"../data/data_cuter/CVXQP2_L","m3","l_bounds","false";
"../data/data_cuter/CVXQP3_L","m3","l_bounds","false";
"../data/data_cuter/DTOC3","m","l_bounds","false";
"../data/data_cuter/HUES-MOD","m","l_eq","true;"
"../data/data_cuter/HUESTIS","m","l_eq","true";
"../data/data_cuter/STCQP1","m","l_bounds","false";
"../data/data_cuter/STCQP2","m","l_bounds","false";
"../data/data_cuter/UBH1","m","l_bounds","false";
"../data/data_cuter/BOYD1","m","l_eq","true"];
end


