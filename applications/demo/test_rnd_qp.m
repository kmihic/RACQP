%


function test_rnd_qp(model_version, n_var, n_blocks, beta, verbose)
  addpath('../../solver/racqp');
  addpath('../../solver/utils');
  addpath('../utils');

  RUNTIME = 120;
  
  disp(" ")
  disp("This version of RACQP implements")
  disp("   min  x'H'x + c'x")
  disp("   s.t. Aeq x = beq")
  disp("        Aineq x <= bineq")
  disp("        lb <= x <= ub")
  disp("        x integer, continuous")
  disp("")

  % model parameters
  if(model_version == 1) %QP
    construct.n_var= n_var;
    construct.Q_sparsity= 0.8;
    construct.kappa = 0;
    construct.c_sparsity= 0;
    construct.Aeq_sparsity= 0.8;
    construct.Aineq_sparsity= 0.8;
    construct.Aineq_n_row= 10;
    construct.Aeq_n_row= 50;
    construct.rnd_seed= 123;
  elseif(model_version == 2) %QP
    construct.n_var= n_var;
    construct.Q_sparsity= 0.1;
    construct.kappa = 0;
    construct.c_sparsity= 0;
    construct.Aeq_sparsity= 0.9;
    construct.Aineq_sparsity= 1;
    construct.Aineq_n_row= 0;
    construct.Aeq_n_row= floor(n_var/2);
    construct.rnd_seed= 123;
  elseif(model_version == 3) %LP
    construct.n_var= n_var;
    construct.Q_sparsity= -1; %if <0, then LP (Q=sparse(n,n))
    construct.kappa = 0;  %regularization kappa*|x|^2
    construct.c_sparsity= 0;
    construct.Aeq_sparsity= 0.8;
    construct.Aineq_sparsity= 0.8;
    construct.Aineq_n_row= 10;
    construct.Aeq_n_row= 50;
    construct.rnd_seed= 123;
  elseif(model_version == 4)  %LP 
    construct.n_var= n_var;
    construct.Q_sparsity= -1;  %if <0, then LP (Q=sparse(n,n))
    construct.kappa = 0; %regularization kappa*|x|^2
    construct.c_sparsity= 0;
    construct.Aeq_sparsity= 0.9;
    construct.Aineq_sparsity= 1;
    construct.Aineq_n_row= 0;
    construct.Aeq_n_row= floor(n_var/2);
    construct.rnd_seed= 123;
  else
    error("Model versions 1, 2 and 3 accepted only")
  end

  problem.construct = construct;
  problem.type = 'construct';

  % RACQP run parameters
  rp.n_blocks = n_blocks; %more blocks (>=2), RP behaves more similar to RAC;
  rp.beta = beta;
  rp.epsilon = 1e-4;
  rp.max_iter = 1500;
  rp.max_rtime = RUNTIME;

  solutions = [];

  solver_mode = "RAC";
  disp("Solving the model using RAC-ADMM")
  [s,model] = demo_rnd_QP(problem,rp,verbose,solver_mode);
  problem.model = model;
  problem.type = 'model';
  s.name = solver_mode;
  solutions = [solutions,s];
  if(n_blocks > 1)
    solver_mode = "RP_ADMM";
    disp("Solving the model using RP-ADMM")
    s = demo_rnd_QP(problem,rp,verbose,solver_mode);
    s.name = solver_mode;
    solutions = [solutions,s];

    solver_mode = "CYCLIC_ADMM";
    disp("Solving the model using cyclic ADMM")
    s = demo_rnd_QP(problem,rp,verbose,solver_mode);
    s.name = solver_mode;
    solutions = [solutions,s];

    msg = "";  
  else   
    msg = "Single block ADMM used. No difference between ADMM types";  
  end
  print_solutions(solutions,false, false, 'ADMM_version' ,msg);

% %PLOT primal residual values at each iteration
% if(verbose > 0)
% figure
% plot(solution_RAC.res_iter, 'b','LineWidth',1)
% hold on
% plot(solution_RP.res_iter, 'r','LineWidth',1)
% plot(solution_ADMM.res_iter, 'g','LineWidth',1)
% legend("RAC","RP-ADMM","CYCLIC-ADMM")
% xlabel("Iterations")
% ylabel("|Ax-b|_{inf}")
% hold off
end
