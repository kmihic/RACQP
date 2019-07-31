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
%addpath('/opt/gurobi/gurobi752/linux64/matlab/');
addpath('../racqp');
addpath('../utils');

  max_time = 5;%600;
  inst = get_instances();
  data_path = "../data/data_qaplib/";
  solutions = [];
  for ii = 1:length(inst)
    filename = data_path+inst(ii,1);
    disp('Solving: '+filename);
    s = verify_QAPLIB_bin(filename+".dat", max_time,solver);
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
inst=[
"lipa80a";
"lipa80b";
"lipa90a";
"lipa90b";
"sko81";
"sko90";
"sko100a";
"sko100b";
"sko100c";
"sko100d";
"sko100e";
"sko100f";
"tai80a";
"tai80b";
"tai100a";
"tai100b";
"tai150b";
"tho40";
"tho150";
"wil50";
"wil100"];
end
