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

%do some simple model verification

function ok = verify_model(model)
  ok = 1;

  %no Inf values
  ok_s = sum(sum(isinf(model.Q))) + sum(isinf(model.c)) ...
    + sum(sum(isinf(model.Aeq))) + sum(isinf(model.beq)) ...
    + sum(sum(isinf(model.Aineq))) + sum(isinf(model.bineq)) ...
    + isinf(model.const) + sum(isinf(model.x0));
  if(ok_s ~= 0)
    ok = 0;
    return;
  end
  % lb/ub check
  lb_i = find(model.lb == Inf);
  ub_i = find(model.ub == -Inf);
  ok_s = sum(lb_i) + sum(ub_i);
  if(ok_s ~= 0)
    ok = 0;
    return;
  end  
end
