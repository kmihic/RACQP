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


function run_maxcut_test()
%addpath('/opt/gurobi/gurobi752/linux64/matlab/');
addpath('../racqp');
addpath('../utils');

max_time = 600;
inst = get_instances();
  solutions = [];
  for ii = 1:length(inst)
    disp('Solving: '+inst(ii,1));
    s = verify_maxcut(inst(ii), max_time);
    %this is a maximization problem, flipping the obj value
    s.sol_obj_val = -s.sol_obj_val;
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
inst=["../data/data_gset/G1";
"../data/data_gset/G6";
"../data/data_gset/G11";
"../data/data_gset/G14";
"../data/data_gset/G18";
"../data/data_gset/G22";
"../data/data_gset/G27";
"../data/data_gset/G32";
"../data/data_gset/G36";
"../data/data_gset/G39";
"../data/data_gset/G43";
"../data/data_gset/G50";
"../data/data_gset/G51";
"../data/data_gset/G55";
"../data/data_gset/G56";
"../data/data_gset/G58";
"../data/data_gset/G60";
"../data/data_gset/G61";
"../data/data_gset/G63";
"../data/data_gset/G67";
"../data/data_gset/G70";
"../data/data_gset/G77";
"../data/data_gset/G81";];
end
