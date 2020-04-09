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
% Note: this is the variant of the original matlab 
%       function (rbf_kernel.m) augmented with 
%       the gamma variable. The original file uses rbf_sigma,
%       which is accessible only for special modes of operation
%       (for more info check matlab documantation on SVM)
%

function kval = mygauss_mat(u,v,rbf_sigma,varargin)

  global g_gamma;

  kval = exp((-g_gamma)*(repmat(sqrt(sum(u.^2,2).^2),1,size(v,1))...
    -2*(u*v')+repmat(sqrt(sum(v.^2,2)'.^2),size(u,1),1)));
end
