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
function term = terminate(run_p, curr_iter, curr_res, time_start, min_iter)
  %we may get lucky and hit primal residual at the very first couple
  %of iterations, but as we do not check for dual residual, we may get a results
  %far away from the optimal. Using min_iter to prevent early termination

  term = false;
  if(curr_iter >= run_p.max_iter ...
     || (curr_iter > run_p.min_iter && curr_res < run_p.epsilon) ...
     || toc(time_start) >= run_p.max_rtime)
     term = true;
  end
end