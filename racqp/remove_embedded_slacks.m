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

function model = remove_embedded_slacks(model)

  sp = [];
  ix=find(diag(model.Q) == 0)';
  %this could be recoded using "any" function
  for ii = ix
    if(nnz(model.Q(ii,:)) == 0 && nnz(model.Q(:,ii)) == 0 && model.c(ii) == 0 && isinf(model.ub(ii)))
      sp = [sp,ii];
    end
  end
  if(length(sp)==0)
    return
  end
  s_ix = [];
  a_ix_rem = [];
  a_ix_add = [];
  for ii = sp
    ix = find(model.Aeq(:,ii) == 1);
    if(length(ix) == 1 && model.Aeq(ix,ii) == 1)
      s_ix = [s_ix,ii];
      a_ix_rem = [a_ix_rem,ix];
      if(isinf(model.lb(ii))) %have unbounded "slack", remove the constraint
        continue;
      end
      a_ix_add = [a_ix_add,ix];
    end
  end

  model.size = model.size - length(s_ix);
  model.Q(:,s_ix) = [];
  model.Q(s_ix,:) = [];
  model.c(s_ix) = [];
  model.x0(s_ix) = [];
  model.lb(s_ix) = [];
  model.ub(s_ix) = [];

  model.Aeq(:,s_ix) = [];
  model.Aineq(:,s_ix) = [];
  Aineq_s = model.Aeq(a_ix_add,:);
  bineq_s = model.beq(a_ix_add);
  model.Aeq(a_ix_rem,:) = [];
  model.beq(a_ix_rem) = [];

  model.Aineq = [model.Aineq;Aineq_s];
  model.bineq = [model.bineq;bineq_s];

end
