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


function x0 = perturb_qap(x, model,n_ch, res_p)
  
  if(res_p == 0)
    x0 = perturb_swap(x, model, n_ch);
  else
    % in a rare case, we may get here
    % make a new feasible starting point
    x=randperm(n_var);
    x0 = zeros(N,1);
    for ii = 1:n_var
      x0((ii-1)*n_var+x(ii)) = 1;
    end
  end

end