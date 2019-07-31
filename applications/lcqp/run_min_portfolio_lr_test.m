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


function run_min_portfolio_lr_test(solver)
addpath('../racqp');
addpath('../utils');


  inst = get_instances();
  solutions = [];
  data_path = "../data/data_markowitz/";
  for ii = 1:length(inst)
    filename = data_path+inst(ii,1);
    disp('Solving: '+filename);
    if(strcmp(solver,'racqp'))
      r_time = 60;
    else
      r_time = inf;
    end
    s = verify_Markowitz_portfolio_selection_data(filename+".mat",true,r_time,solver);
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
inst=["lowrank_daily";  
"lowrank_monthly";  
"lowrank_quarterly"];
end
