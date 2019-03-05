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

%TEST: QAP relaxed (x\in[0,1])
%      Using QAPLIB benchmark problems 
%      (http://anjos.mgi.polymtl.ca/qaplib/)
%
% Output: qapcont_out


clear all;
%load and convert the problem to rac style
datapath="../data/data_qaplib/";
filename=datapath+"wil50.dat";
%no grouping
%model = get_QAPLIB(filename);
%group vars
model = get_QAPLIB(filename,true);
%This is the relaxed QAP, keep x within bounds, but remove binary constraint
model.binary=[];

%load default runtime parameters 
default_run_params
r = floor(sqrt(model.size));
run_params.n_blocks = r;
run_params.beta = r;
run_params.epsilon = 1e-4;
run_params.epsilon_d = 1e-2;
run_params.max_iter = 1000;

%turn on verbose
run_params.debug = 1;

%run the solver
disp('Running RAC MIP')
qapcont_out = RACQP(model, run_params);

disp('Solution, relaxed QAP:')
qapcont_out






