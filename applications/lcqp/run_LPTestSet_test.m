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


function run_LPTestSet_test(solver)
%addpath('/opt/gurobi/gurobi752/linux64/matlab/');
addpath('../racqp');
addpath('../utils');

  inst = get_instances();
  solutions = [];
  data_path = "../data/data_mps/";
  tmp_dir = "./tmp_rac/";
  mkdir(tmp_dir,'s');
  for ii = 1:length(inst)
    filename = inst(ii,1)+".mps";
    disp('Solving: '+data_path+filename);
    disp("Unzipping the file")
    gunzip(data_path+filename+".gz", tmp_dir);
    s = verify_LPTestSet(tmp_dir+filename,str2num(inst(ii,2)),str2num(inst(ii,3)),inst(ii,4),solver);
    solutions = [solutions,s];
    delete(tmp_dir+filename);
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
  rmdir(tmp_dir,'s')
  T = table(name,obj_val,run_time);
  T.Properties.VariableNames={'Data_file','obj_val', 'run_time'};
  disp(T);

%quit
end

function inst = get_instances()
inst=["wide15", "10", "200", "none";
"square15", "10", "200", "none";
"long15", "10", "200", "none";
"nug30", "200", "200", "l_bounds"];
end


