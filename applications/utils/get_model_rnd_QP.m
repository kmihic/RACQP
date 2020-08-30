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



function model=get_rnd_QP(p_setup, add_bounds, Q_setup)
  %set problem generator seed
  rng(p_setup.rnd_seed);

  %get Q
  %in a controlled way: n_var,cond_num,eta,zeta,sparsity_Q,rank
  if(p_setup.Q_sparsity < 0)
    Q = sparse(p_setup.n_var,p_setup.n_var);
  elseif(nargin > 2 )
    Q = makeQ(p_setup, Q_setup);
  else
    Q = makeQ(p_setup);
  end
  Q = Q + p_setup.kappa * speye(p_setup.n_var,p_setup.n_var);
  %get c
  c = get_c(p_setup.n_var,p_setup.c_sparsity);
  
  
  %get eq constraints
  x0 = rand(length(c),1);
  if(p_setup.Aeq_sparsity==1 || p_setup.Aeq_n_row == 0)
    Aeq=sparse(0,p_setup.n_var);
    beq=sparse(0,1);
  else
    Aeq = get_A(p_setup.n_var,p_setup.Aeq_n_row,p_setup.Aeq_sparsity);  
    beq = sparse(Aeq*x0);
  end
  %get ineq constraints
  if(p_setup.Aineq_sparsity==1 || p_setup.Aineq_n_row == 0)
    Aineq=sparse(0,p_setup.n_var);
    bineq=sparse(0,1);
  else
    Aineq = get_A(p_setup.n_var,p_setup.Aineq_n_row,p_setup.Aineq_sparsity);  
    bineq = sparse(Aineq*x0+ones(p_setup.Aineq_n_row,1));
  end
  
  model.size = p_setup.n_var;
  model.Q = Q;
  model.c = c;
  model.Aeq = Aeq;
  model.beq = beq;  
  model.Aineq = Aineq;
  model.bineq = bineq;
  model.x0 = zeros(p_setup.n_var,1);
  model.integers = [];
  model.binary = [];
  if(nargin > 1 && add_bounds)
    %x0 is from U(0,1)
    model.lb = zeros(p_setup.n_var,1); 
    model.ub = ones(p_setup.n_var,1);
  else
    model.lb = ones(p_setup.n_var,1)*(-Inf);  
    model.ub = ones(p_setup.n_var,1)*(Inf);
  end  
  model.const = 0; 
end

function c=get_c(n_var,sparsity)
  if(sparsity == 1)
    c = sparse(n_var,1);
    return;
  end
  max_c=0;
  c=zeros(n_var,1);

  for j=1:1:n_var
    if rand<sparsity
      continue;
    end
    %val=normrnd(0,1);
    val = rand();
    c(j)=val;
    max_c=max(max_c,val);
  end
  c=sparse(c./max_c);
  
end

function A=get_A(n_var,n_row,sparsity)
  max_A=0;
  A=zeros(n_row,n_var);
  for i=1:1:n_row
    for j=1:1:n_var
      if rand<sparsity
        continue;
      end
      val=normrnd(0,1);
      A(i,j)=val;
      max_A=max(max_A,val);
    end
  end

  A=sparse(A./max_A);  
end


function Q = makeQ(p_setup,Q_setup)

  if(nargin == 2)
    Q = get_controlled_Q(p_setup.n_var,p_setup.Q_sparsity,Q_setup);
  else
    Q = get_rnd_Q(p_setup.n_var, p_setup.Q_sparsity);
  end
  Q=sparse(Q);
end

% Note: sparsity must be very high for sprand to produce correct
%       results. See help sprand
function Q=get_rnd_Q(n_var, sparsity)
  if(sparsity==1)
    Q=sparse(n_var,n_var);
    return;
  end
  A=sprand(n_var,n_var,1-sparsity);
  Q = A+A';
  %Q has all positive values (rand>=0), no need for abs
  %make it psd
  off_diag = sum(Q,2) - diag(Q) +1e-10; 
  max_off_diag = max(off_diag);
  Q = Q+diag(off_diag);
  %normalize Q
  Q = Q./max_off_diag;
end

function Q=get_controlled_Q(n_var,sparsity,s)
  if(sparsity==0)
    Q=sparse(n_var,n_var);
    return;
  end
  cond_num = s.cond_num;
  eta = s.eta;
  zeta = s.zeta;
  %%zeta>0, eta>0
  if (zeta <= 0)
   error('Error. Input zeta must be > 0 ')
  end
  if (eta < 0 || eta > 1)
   error('Error. Input eta must be [0,1] ')
  end
  %%Generate diagonal matrix 
  S=exp(rand(1,n_var)); 
  S=S./max(S);

  %cond_num: 
  [S_min,S_min_index]=min(S);
  S(S_min_index)=1/cond_num;

  %%generate V, Vt
  V=eta*rand(n_var,n_var);
  V=V+(1-eta)*eye(n_var);
  Vt=V;
  V=V.*repmat(S,n_var,1);

  Q=sparse(n_var,n_var);
  max_q=0;
  sum_row=zeros(n_var,1);
  for i=1:1:n_var
      for j=i:1:n_var
          if (i~=j) && (rand < sparsity)
              continue;
          end
          val=zeta+V(i,:)*transpose(Vt(j,:));
          Q(i,j)=val;
          max_q=max(abs(val),max_q);
          if (sparsity>0) && (i~=j)
              sum_row(i)=sum_row(i)+val;
              sum_row(j)=sum_row(j)+val;
          end
      end
      if (sparsity>0) && (Q(i,i)<sum_row(i))
           Q(i,i)=Q(i,i)+sum_row(i);
           max_q=max(max_q,abs(Q(i,i)));
      end 
  end
  Q=Q./max_q;
  Q=Q+Q'-diag(diag(Q));

  Q=sparse(Q);
end
