

function x0 = get_init_pt(lb, ub)

  n = length(lb);
  x0 = zeros(n,1);
  lb_ix = find(lb>0);
  ub_ix = find(ub<0);
  x0(lb_ix) = lb(lb_ix);
  x0(ub_ix) = ub(ub_ix);
end
