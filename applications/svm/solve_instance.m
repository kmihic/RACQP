
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
function sol = solve_instance(model_libsvm, solver, rac_rp, no_test,beta)
  
  sol = [];
  if(strcmpi(solver,'racqp'))
    for b = beta
      rac_rp.beta = b;
      s = use_rac(model_libsvm,rac_rp, no_test);
      if(length(beta) > 1)
        s.beta = b;
      end
      s.C = model_libsvm.C;
      s.sigma = model_libsvm.sigma;
      sol=[sol s];
    end
  elseif(strcmpi(solver,'libsvm'))
    %run libsvm    
    sol = use_libsvm(model_libsvm, no_test);
    sol.C = model_libsvm.C;
    sol.sigma = model_libsvm.sigma;
  elseif(strcmpi(solver,'matlab'))
    %run matlab
    sol = use_matlab(model_libsvm, no_test);
    sol.C = model_libsvm.C;
    sol.sigma = model_libsvm.sigma;
  end
    
end
