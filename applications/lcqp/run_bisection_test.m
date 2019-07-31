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


function run_maxcut_test(r_time,solver)
%addpath('/opt/gurobi/gurobi752/linux64/matlab/');
addpath('../racqp');
addpath('../utils');

  max_time = r_time;
  inst = get_instances();
  solutions = [];  
  data_path = "../data/data_gset/";
  for ii = 1:length(inst)
    filename = data_path+inst(ii,1);
    disp('Solving: '+filename);
    s = verify_bisection(filename, max_time, solver);
    %this is a maximization problem, flipping the obj value
    s.sol_obj_val = -s.sol_obj_val;
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
end

function inst = get_instances()
inst=["G1";
"G6";
"G11";
"G14";
"G18";
"G22";
"G27";
"G32";
"G36";
"G39";
"G43";
"G50";
"G51";
"G55";
"G56";
"G58";
"G60";
"G61";
"G63";
"G67";
"G70";
"G77";
"G81"];
end
