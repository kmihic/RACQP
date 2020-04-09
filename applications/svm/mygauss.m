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

function K_rbf = mygauss(u,v)

  global g_gamma;

  if(nargin == 2)
    u_sq=sum(u.^2,2);
    v_sq=sum(v.^2,2);
    K=bsxfun(@minus,u_sq,(2*u)*v.');
    K=bsxfun(@plus,v_sq.',K);
    K_rbf=exp((-g_gamma)*K);
  else
    utu = u*u';
    du = diag(utu);
    K_rbf=exp(-g_gamma*(du-2*utu+du'));
  end
end
