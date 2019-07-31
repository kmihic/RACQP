
function label = classify(u,v, y_alpha, b)

  %saving memory at the expense of computational time
  block_size = 1000;
  N = size(u,1);
  n_blocks = ceil(N/block_size);
  label = [];
  %u_sq=sum(u.^2,2);
  %v_sq=sum(v.^2,2);
  for jj = 0:(n_blocks-1)
    bb = jj*block_size;
    be = min(bb+block_size, N);
    ii = (bb+1):be;
    K_rbf = mygauss(u(ii,:),v);
    label(ii) = sign(K_rbf*y_alpha + b);
  end
  label = label';
end