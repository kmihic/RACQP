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

function model = clean_model_constraints(model, remove_slacks)

  if(nargin<=1)
   remove_slacks = false;
  end
  
  lb_inf = isinf(model.lb);
  ub_inf = isinf(model.ub);
  no_bounds = lb_inf & ub_inf;

  no_coeff_c = model.c == 0;
  if(~issymmetric(model.Q))
    no_coeff_Q_cols = sum(model.Q,1)==0;
    no_coeff_Q_rows = sum(model.Q,2)==0;
    no_coeff_Q = no_coeff_Q_cols' & no_coeff_Q_rows;
  else
    no_coeff_Q = (sum(model.Q,1)==0)';
  end
  
  % find unbounded variables used in constraints only once -> those rows
  % are always feasible if vars are not operands of f_0 (i.e coeffs are zero)
  no_f0 = no_coeff_c & no_coeff_Q;
  no_bc = no_bounds & no_f0; 
  no_bc_x_ix = find(no_bc);
  row_no_bc_eq = model.Aeq(:,no_bc_x_ix)~=0;
  row_no_bc_ineq = model.Aineq(:,no_bc_x_ix)~=0;
  single_use = find(sum([row_no_bc_eq;row_no_bc_ineq],1)==1);
  rm_x_ix = no_bc_x_ix(single_use);
  [invariant_rows_eq_ix, ~] = find(row_no_bc_eq(:,single_use));
  [invariant_rows_ineq_ix, ~] = find(row_no_bc_ineq(:,single_use));
  % in a case we have multiple cols 
  rm_rows_eq = unique(invariant_rows_eq_ix);
  rm_rows_ineq = unique(invariant_rows_ineq_ix);

  % find embedded slacks and re-create ineqality constraints
  if(remove_slacks)
    x_ge_zero = (model.lb == 0) & ub_inf & no_f0;
    x_le_zero = (model.ub == 0) & lb_inf & no_f0;
    x_ge_zero_ix = find(x_ge_zero);
    x_le_zero_ix = find(x_le_zero);
    row_ge_zero = model.Aeq(:,x_ge_zero_ix)~=0;    
    row_le_zero = model.Aeq(:,x_le_zero_ix)~=0;
    single_use_ge = find(sum(row_ge_zero,1)==1);
    single_use_le = find(sum(row_le_zero,1)==1);
    rm_x_ix = [rm_x_ix;x_ge_zero_ix(single_use_ge);x_le_zero_ix(single_use_le)];

    [slack_rows_ge_pos, ~] = find(row_ge_zero(:,single_use_ge)>0); % Ax +s = b, s >=0
    [slack_rows_ge_neg, ~] = find(row_ge_zero(:,single_use_ge)<0); % Ax -s = b, s >=0
    [slack_rows_le_pos, ~] = find(row_le_zero(:,single_use_le)>0); % Ax +s = b, s <=0
    [slack_rows_le_neg, ~] = find(row_le_zero(:,single_use_le)<0); % Ax -s = b, s <=0
    slack_rows_ge_pos = unique(slack_rows_ge_pos);
    slack_rows_ge_neg = unique(slack_rows_ge_neg);
    slack_rows_le_pos = unique(slack_rows_le_pos);
    slack_rows_le_neg = unique(slack_rows_le_neg);

    % Ax <= b: Ax+s=b, s>= 0; Ax-s=b, s<=0
    % Ax >= b: Ax-s=b, s >= 0; Ax+s=b, s <= 0    
    add_rows_ineq = [slack_rows_ge_pos;slack_rows_ge_neg;
                     slack_rows_le_pos;slack_rows_le_neg];
    neg_ix = [slack_rows_le_pos;slack_rows_ge_neg];
    model.Aeq(:,neg_ix)=-model.Aeq(:,neg_ix);
    model.beq(neg_ix)=-model.beq(neg_ix);
  
    % move inequality rows from Aeq to Aineq
    model.Aineq = [model.Aineq; model.Aeq(add_rows_ineq,:)];
    model.bineq = [model.bineq; model.beq(add_rows_ineq)];
    rm_rows_eq = [rm_rows_eq;add_rows_ineq];
  end
  
  % remove the rows
  model.Aeq(rm_rows_eq,:) = [];
  model.beq(rm_rows_eq) = [];  
  model.Aineq(rm_rows_ineq,:) = [];
  model.bineq(rm_rows_ineq) = [];
  % remove the columns
  model.size = model.size - length(rm_x_ix);
  model.Q(:,rm_x_ix) = [];
  model.Q(rm_x_ix,:) = [];
  model.c(rm_x_ix) = [];
  model.x0(rm_x_ix) = [];
  model.lb(rm_x_ix) = [];
  model.ub(rm_x_ix) = [];
  model.Aeq(:,rm_x_ix) = [];
  model.Aineq(:,rm_x_ix) = [];


end
