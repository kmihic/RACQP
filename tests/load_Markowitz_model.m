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


function model = load_Markowitz_model(filename, p_cut, kappa, low_rank, mip)

load(filename);
%get the problem 
binary=[];
if(low_rank)
  c = c_lowrank;
  N=length(c);  
  Q = kappa*speye(N);
  Q(index_continuous,index_continuous)=1;
  beq = b_lowrank;
  Aeq = A_lowrank;  
  %x vars are either binary or boxed cont
  lb = zeros(N,1);
  lb(index_continuous) = -inf;    
  ub = ones(N,1);
  ub(index_continuous) = inf;  
  if(mip)
    n_x = length(index_integer);
    beq(1) = min(n_x,ceil(p_cut*n_x));    
    binary = index_integer;        
  else
    beq(1) = 1;
  end
else
  N=size(Q,1);
  Q=sparse(Q) + kappa*speye(N);
  %x vars are either binary or boxed cont
  lb = zeros(N,1);  
  ub = ones(N,1);
  Aeq = sparse(ones(1,N)); 
  if(mip)
    beq = min(N,ceil(p_cut*N));
    binary = (1:N)';    
  else
    beq = 1;    
  end
end


model.size = N;
model.c = c;
model.Q = Q;
model.Aeq = Aeq;
model.beq = beq;  
model.Aineq = sparse(0,N);
model.bineq = [];
model.x0 = zeros(N,1);
model.integers = [];
model.binary = binary;
model.lb = lb;  
model.ub = ub;
model.const = 0;

end
