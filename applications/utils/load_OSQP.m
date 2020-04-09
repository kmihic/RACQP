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

function model = load_OSQP(filename)
  %loads OSQP model and converts it to RACQP format

  load(filename);
  osqp.A = A;
  osqp.P = P;
  osqp.n = n;
  osqp.m = m;
  osqp.q = q;
  osqp.u = u;
  osqp.l = l;
  osqp.r = r;
  model = osqp_to_rac(osqp);
end

function model = osqp_to_rac(osqp)

  n = osqp.n;
  model.size = n;
  model.Q = osqp.P/2;
  model.c = sparse(osqp.q);

  Aeq=sparse(0,n);
  beq=sparse(0,1);
  Aineq_le=sparse(0,n);
  Aineq_ge=sparse(0,n);
  bineq_le=sparse(0,1);
  bineq_ge=sparse(0,1);

  %make sure every var has a bound, +-inf (1e10), will correct later
  U=ones(n,1)*Inf;
  L=-ones(n,1)*Inf;
  %correct osqp lhs/rhs
  ix = find(osqp.u==1e20);
  osqp.u(ix) = Inf;
  ix = find(osqp.l==-1e20);
  osqp.l(ix) = -Inf;
  ae_ix = sparse(size(osqp.A,1),1);
  ail_ix = sparse(size(osqp.A,1),1);
  aig_ix = sparse(size(osqp.A,1),1);
  for ii=1:size(osqp.A,1)
    ix= find(osqp.A(ii,:) ~= 0); 
    if(length(ix) == 1) %have one non-zero component, thus bounds are defined here
      if(~isinf(osqp.u(ii)))
        U(ix) = osqp.u(ii)/osqp.A(ii,ix);
      end
      if(~isinf(osqp.l(ii)))
        L(ix) = osqp.l(ii)/osqp.A(ii,ix);
      end
    else  %we have a constraint
      if(isinf(osqp.l(ii)) && isinf(osqp.u(ii)))
        continue; %drop it; error in the model
      elseif(osqp.l(ii)==osqp.u(ii)) %equality 
        ae_ix(ii) = 1;
        %Aeq = sparse([Aeq;osqp.A(ii,:)]);
        %beq = [beq;osqp.u(ii)];
      else  %inequality 
        if(~isinf(osqp.u(ii)))
          ail_ix(ii) = 1;
          %Aineq_le = [Aineq_le;osqp.A(ii,:)];
          %bineq_le = [bineq_le;osqp.u(ii)];
        end
        if(~isinf(osqp.l(ii)))
          aig_ix(ii) = 1;
          %Aineq_ge = [Aineq_ge;osqp.A(ii,:)];
          %bineq_ge = [bineq_ge;osqp.l(ii)];
        end
      end
    end
  end
  ix = find(ae_ix);
  beq = osqp.u(ix);
  Aeq = osqp.A(ix,:);
  ix = find(ail_ix);
  bineq_le = osqp.u(ix);
  Aineq_le = osqp.A(ix,:);
  ix = find(aig_ix);
  bineq_ge = osqp.l(ix);
  Aineq_ge = osqp.A(ix,:);
  
  model.Aeq = Aeq;
  model.beq = sparse(beq);
  model.Aineq = [Aineq_le;-Aineq_ge];
  model.bineq = sparse([bineq_le;-bineq_ge]);
  model.x0 = get_init_pt(model.lb, model.ub);%get_init_x(L,U);
  model.integers = [];
  model.binary = [];
  model.lb = L;
  model.ub = U;
  model.const = osqp.r;

end

%find init point at rnd, within bouds
function x = get_init_x(L,U)

  rng(1);
  x=zeros(length(L),1);
  for ii = 1:length(L)
    if( isinf(L(ii)) && isinf(U(ii)) ) %no limits
      x(ii) = randn(1);
    elseif( isinf(U(ii)) )  %have a lower bound
      x(ii) = L(ii)+abs(randn(1));
    elseif( isinf(L(ii)) ) %have an upper bound
      x(ii) = U(ii) - abs(randn(1));
    else %have both bounds
      x(ii) = rand(1)*(U(ii)-L(ii))+L(ii);
    end
  end
end

  
