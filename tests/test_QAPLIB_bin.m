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

%TEST: QAP binary (x\in{0,1})
%      Using QAPLIB benchmark problems 
%      (http://anjos.mgi.polymtl.ca/qaplib/)
%
% Output: qapbin_out

clear all
%load and convert the problem to rac style
datapath="../data/data_qaplib/";
filename=datapath+"wil50.dat";
%group vars
%fileaname="tai10a.dat";
model = get_QAPLIB(filename,true);

nsize = model.size;
lambda = 0.3;
max_time = 30;
p_trial = 0.001;


%load default runtime parameters 
default_run_params
r = floor(sqrt(model.size));
run_params.n_blocks = r;
run_params.beta = model.size;
run_params.epsilon = 1e-4;
run_params.epsilon_d = -Inf;
run_params.max_iter = 1000;


%turn off verbose
run_params.debug = 0;

%load default gurobi parameters
default_gurobi_parameters
%change how long gurobi will spend on each subproblem
gurobi_params.TimeLimit = 5;
 
%gurobi_params.outputflag=1;
run_params.gurobi_params = gurobi_params;

%set MIP run parameters

run_mip.permute_type = 'user_defined';
run_mip.permute_f = @perturb_qap;
run_mip.permute_dist = 'exponential'; 
run_mip.permute_min = 0.005*nsize; 
run_mip.permute_max = 0.7*nsize; 
run_mip.permute_mu = lambda*nsize; 
run_mip.max_iter = 10000000; %time limited 
run_mip.max_nperturb = 10000000; %time limited
run_mip.max_rtime = max_time; 
run_mip.rnd_seed = 123; 
run_mip.n_perturb_trial = max(2,round(p_trial*nsize));
run_mip.run_sub = run_params;
run_mip.debug = 1;

%print run parameters
disp("RUN PARAMETERS")
disp("run_mip")
disp(run_mip)
disp("run_sub")
disp(run_mip.run_sub)
disp("gurobi_parameters")
disp(run_mip.run_sub.gurobi_params)

%run the solver
disp('Running RAC MIP')
qapbin_out = RACQP(model, run_mip);

disp('Solution, binary QAP:')
qapbin_out






