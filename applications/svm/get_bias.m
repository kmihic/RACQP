
function b = get_bias(u,y_alpha, y_label)

  %saving memory at the expense of computational time
  block_size = 1000;
  N = length(y_alpha);
  n_blocks = ceil(N/block_size);
  b = [];
  for jj = 0:(n_blocks-1)
    bb = jj*block_size;
    be = min(bb+block_size, N);
    ii = (bb+1):be;
    K_rbf = mygauss(u,u(ii,:));
    b(ii) = y_label(ii) - (y_alpha'*K_rbf)';
  end
  b = mean(b);
end