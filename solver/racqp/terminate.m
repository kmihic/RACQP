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
function term = terminate(run_p, curr_iter, curr_res, time_start, mip)
  term = false;
  if(curr_iter >= run_p.max_iter ...
     || (curr_iter > run_p.min_iter ...
           && curr_res.prim < run_p.epsilon_prim ...
           && curr_res.dual < run_p.epsilon_dual) ...
     || (mip && curr_iter > run_p.min_iter ...
           && curr_res.prim < run_p.epsilon_prim) ...
     || toc(time_start) >= run_p.max_rtime)
     term = true;
  end
end
