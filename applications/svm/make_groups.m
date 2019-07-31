function groups = make_groups(n_blocks, N)

  block_size = floor(N/n_blocks);
  o_blocks = N-block_size*n_blocks;
  b_size = ones(n_blocks,1) * block_size;
  b_size(1:o_blocks) = b_size(1:o_blocks)+1;

  x_ix_perm = 1:N;%randperm(N);
  groups = cell(n_blocks,2);
  bb = 0;
  for block_ix = 1:n_blocks
    block_size = b_size(block_ix);
    be = bb+block_size;
    bnew= x_ix_perm((bb+1):be);
    groups{block_ix,1} = bnew;
    bb = be;
  end
 

end