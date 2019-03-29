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


% Verify results: QAP binary (x\in{0,1})
%      Using QAPLIB benchmark problems 
%      (http://anjos.mgi.polymtl.ca/qaplib/)
%
% Output: qapbin_out


function [qapbin_out, model] = verify_QAPLIB_bin(filename, max_time) 
%addpath('/opt/gurobi/gurobi752/linux64/matlab/');
addpath('../racqp');
addpath('../utils');

%hardcode experimental run-parameters
lambda = 0.4;
p_trial = 0.005;

%load and convert the problem to rac style
disp("LOADING THE MODEL...")
model = load_QAPLIB(filename);
nsize = model.size;
%do not split groups
model.group_mode = 'RAC_NO_SPLIT';

%load default runtime parameters and change them
default_run_params
r = floor(sqrt(model.size));
run_params.n_blocks = ceil(r/2);
run_params.beta = model.size;
run_params.max_iter = 1000;
%QAP are dense and sparse problems, but Gurobi requires sparse Hessian, and
%we use Gurobi to handle binary problems. No benefit
%of changing full->sparse -> full,...

%turn off verbose
run_params.debug = 0;

%load default gurobi parameters
default_gurobi_parameters
%change max time gurobi can spend on each subproblem
gurobi_params.TimeLimit = 5;
 
%gurobi_params.outputflag=1;
run_params.gurobi_params = gurobi_params;

%set MIP run parameters

run_mip.permute_type = 'user_defined';
run_mip.permute_f = @perturb_qap;
run_mip.permute_dist = 'exponential'; 
run_mip.permute_min = 2; 
run_mip.permute_max = r; 
run_mip.permute_mu = lambda*r; 
run_mip.max_iter = 10000000; %time limited 
run_mip.max_nperturb = 10000000; %time limited
run_mip.max_rtime = max_time; 
run_mip.rnd_seed = 123; 
run_mip.n_perturb_trial = max(2,round(p_trial*nsize));
run_mip.run_sub = run_params;
run_mip.debug = 1;
run_mip.mip_epsilon = 0;

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
disp(qapbin_out)

%quit
end






