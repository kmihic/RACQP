function kval = mygauss_mat(u,v,rbf_sigma,varargin)
%RBF_KERNEL Radial basis function kernel for SVM functions

% Copyright 2004-2012 The MathWorks, Inc.

global g_gamma;
 
if nargin < 3 || isempty(rbf_sigma)
    rbf_sigma = 1;
else
     if ~isscalar(rbf_sigma) || ~isnumeric(rbf_sigma)
        error(message('stats:rbf_kernel:RBFSigmaNotScalar'));
    end
    if rbf_sigma == 0
        error(message('stats:rbf_kernel:SigmaZero'));
    end
    
end

kval = exp((-g_gamma)*(repmat(sqrt(sum(u.^2,2).^2),1,size(v,1))...
    -2*(u*v')+repmat(sqrt(sum(v.^2,2)'.^2),size(u,1),1)));
