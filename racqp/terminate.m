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
function term = terminate(run_p, curr_iter, curr_res, time_start)

  term = false;
  if(curr_iter >= run_p.max_iter ...
     || (curr_res < run_p.epsilon) ...
     || toc(time_start) >= run_p.max_rtime)
     term = true;
  end
end