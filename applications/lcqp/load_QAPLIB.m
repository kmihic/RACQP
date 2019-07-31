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


function model = load_QAPLIB(filename)
  delimiter = ' ';
  headLineN = 1;
  
  fi=fopen(filename,'r');
  data_in = fscanf(fi,'%d');
  fclose(fi);
  n_var = data_in(1);
  N = n_var^2;
  F = reshape(data_in(2:N+1),n_var,[]);
  D = reshape(data_in(N+2:2*N+1),n_var,[]);

  if(length(data_in) ~= 2*N+1)
    error('Error. Something is wrong with input data');
  end  


  Q = kron(D,F);
  %normalize Q
  max_Q = max(abs(Q(:)));
  Q = Q./max_Q;
  %check if symmetric
  %some instances may be asymmetric (e.g. lipa*)
  %if(~issymmetric(Q))  
  %  error('Error. Something is wrong with input data: matrices not symmetric')
  %end
  %need Q to be PSD, making it diagonally dominant 
  Q_diag = 2 * max( (sum(abs(Q),2)-diag(abs(Q))) );
  Q = Q+Q_diag*eye(N);
  
  %make constraints
  %create A to check for doubly stochastic matrix X
  Aeq = zeros(2*n_var,N);
  for ii = 1:n_var
    Aeq(ii,(ii-1)*n_var+1:ii*n_var) = ones(1,n_var);
  end
  
  for ii = 1:n_var
    Aeq(n_var+ii,[ii:n_var:(n_var-1)*n_var+ii]) = ones(1,n_var);
  end
  beq = ones(2*n_var,1);
  % make a feasible starting point
  x=randperm(n_var);
  x0 = zeros(N,1);
  for ii = 1:n_var
    x0((ii-1)*n_var+x(ii)) = 1;
  end
 
  model.size = N;
  model.Q = sparse(Q);
  model.c = sparse(N,1);

  %model.Aeq = sparse(Aeq(1+n_var:2*n_var,:));%sparse(Aeq);
  %model.beq = sparse(beq(1+n_var:2*n_var));  

  model.Aeq = sparse(Aeq);
  model.beq = sparse(beq);
  model.Aineq = sparse(0,N);
  model.bineq = sparse(0,1);
  model.x0 = x0;
  model.integers = [];
  model.binary = (1:N)';
  model.lb = zeros(N,1);  
  model.ub = ones(N,1);
  model.const = -Q_diag*n_var;

  %model.local_constraints.Aeq = sparse(Aeq(1:n_var,:));%sparse(0,N);
  %model.local_constraints.beq = sparse(beq(1:n_var));%sparse(0,1);

  model.local_constraints.Aeq = sparse(0,N);
  model.local_constraints.beq = sparse(0,1);
  model.local_constraints.Aineq = sparse(0,N);
  model.local_constraints.bineq = sparse(0,1);
  model.local_constraints.lb=zeros(N,1);  
  model.local_constraints.ub=ones(N,1);  
   
  model.norm_coeff = max_Q;

  model.groups = group_vars(Aeq);  
    
end

% For QAP blocks can be either constructed considering first n
% constraints that test for sum row = 1 of using the second n constraints
% enforcing sum col = 1. 
% Using the first n constraints. No overlaps.
function groups = group_vars(A)
  
  groups={};
  for ii = 1:size(A,1)/2
    ix = find(A(ii,:));
    groups(ii,:) = {ix, []};
  end   

end



