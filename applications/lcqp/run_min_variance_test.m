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


function run_min_variance_test(solver)
addpath('../racqp');
addpath('../utils');

  inst = get_instances();
  solutions = [];
  data_path = '../data/data_markowitz/';
  for ii = 1:length(inst)
    filename = data_path+inst(ii,1);
    disp('Solving: '+filename);
    s = verify_Markowitz_min_variance_data(filename+".mat",inst(ii,2),solver);
    % low rank model uses single block, localized kkt, no dual var
    % avaiable for dual_res. Adding it here so matlab will not complain for 
    % different structures
    s.sol_residue_d = nan;
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
inst=["regular_quarterly","false";
"regular_monthly","false";
"regular_daily","false"; 
"lowrank_quarterly","true";
"lowrank_monthly","true";  
"lowrank_daily","true"];
end
