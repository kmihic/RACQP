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


function run_qap_relaxed_test()
addpath('../racqp');
addpath('../utils');

inst = get_instances();
  solutions = [];
  for ii = 1:length(inst)
    disp('Solving: '+inst(ii,1));
    s = verify_QAPLIB_cont(inst(ii));
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
inst=["../data/data_qaplib/dre110.dat"; 
"../data/data_qaplib/esc128.dat"; 
"../data/data_qaplib/dre132.dat"; 
"../data/data_qaplib/sko100a.dat"; 
"../data/data_qaplib/sko100b.dat"; 
"../data/data_qaplib/sko100c.dat"; 
"../data/data_qaplib/sko100d.dat"; 
"../data/data_qaplib/sko100e.dat"; 
"../data/data_qaplib/sko100f.dat"; 
"../data/data_qaplib/tai100a.dat"; 
"../data/data_qaplib/tai100b.dat"; 
"../data/data_qaplib/tai125e01.dat"; 
"../data/data_qaplib/tai150b.dat"; 
"../data/data_qaplib/tho150.dat"; 
"../data/data_qaplib/wil100.dat"];
end
