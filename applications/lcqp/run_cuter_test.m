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


function run_cuter_test(solver)
%addpath('/opt/gurobi/gurobi752/linux64/matlab/');
addpath('../racqp');
addpath('../utils');

  inst = get_instances();
  solutions = [];
  data_path = "../data/data_cuter/";
  for ii = 1:length(inst)
    filename = data_path+inst(ii,1);
    disp('Solving: '+filename);
    s = verify_CUTEr(filename,inst(ii,2),inst(ii,3),inst(ii,4), solver);
    solutions = [solutions,s];
  end
  disp(" ")
  disp("#####################")
  disp('SUMMARY')
  name = [];
  obj_val = [];
  run_time = [];
  for ii = 1: length(solutions)
    name = [name;inst(ii,1)];
    obj_val = [obj_val;solutions(ii).sol_obj_val];
    run_time = [run_time;solutions(ii).rac_time];
  end
  T = table(name,obj_val,run_time);
  T.Properties.VariableNames={'Data_file','obj_val', 'run_time'};
  disp(T); 
  
%quit
end


function inst = get_instances()
inst=["AUG2D","m","none","false";
"AUG2DC","m","none","false";
"AUG2DCQP","m","l_bounds","false";
"AUG2DQP","m","l_bounds","false";
"AUG3D","m","none","false";
"AUG3DC","m","none","false";
"AUG3DCQP","m","l_bounds","false";
"AUG3DQP","m","l_bounds","false";
"BOYD1","m","l_eq","true";
"CONT-050","m","l_bounds","false";
"CONT-100","m","l_bounds","false";
"CONT-101","m","l_bounds","false";
"CONT-200","m","l_bounds","false";
"CONT-201","m","l_bounds","false";
"CONT-300","m","l_bounds","false";
"CVXQP1_L","m3","l_bounds","false";
"CVXQP2_L","m3","l_bounds","false";
"CVXQP3_L","m3","l_bounds","false";
"DTOC3","m","l_bounds","false";
"HUES-MOD","m","l_eq","true;"
"HUESTIS","m","l_eq","true";
"STCQP1","m","l_bounds","false";
"STCQP2","m","l_bounds","false";
"UBH1","m","l_bounds","false"];
end


