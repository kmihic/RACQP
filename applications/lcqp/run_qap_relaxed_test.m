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


function run_qap_relaxed_test(solver)
addpath('../racqp');
addpath('../utils');

  inst = get_instances();
  solutions = [];
  data_path = '../data/data_qaplib/';
  for ii = 1:length(inst)
    filename = data_path+inst(ii);
    disp('Solving: '+filename);
    s = verify_QAPLIB_cont(filename+".dat", solver);
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
inst=["dre110"; 
"esc128"; 
"dre132"; 
"sko100a"; 
"sko100b"; 
"sko100c"; 
"sko100d"; 
"sko100e"; 
"sko100f"; 
"tai100a"; 
"tai100b"; 
"tai125e01"; 
"tai150b"; 
"tho150"; 
"wil100"];
end
