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

function x0 = perturb_swap(x, model, n_ch)
  x0 = x;
  x_green = find(~x);
  x_red = find(x);
  
  %chech how many we can swap
  n_ch = min(n_ch, min(length(x_green), length(x_red)));
  %find vars to swap
  s_green = datasample(x_green,n_ch,'Replace',false);  
  s_red = datasample(x_red,n_ch,'Replace',false);
  %swap the values
  x0(s_green)=x(s_red);
  x0(s_red)=x(s_green);
end