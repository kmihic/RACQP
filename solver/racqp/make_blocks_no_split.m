function out = make_blocks_no_split(n_vars,groups, block_sizes, no_split)
  
  n_group = size(groups,1);
  var_used = zeros(n_vars,1);
  n_blocks = length(block_sizes);
  blocks = cell(1,1);
  grp_order = randperm(n_group);
  cnt_b = 1;
  block = [];
  n_sum = 0;

  space_left = block_sizes(cnt_b);
  while(n_group > 0)
    g_ix = grp_order(n_group); 
    v_ix = groups{g_ix,1};
    % remove already assigned vars
    v_ix = v_ix(find(~var_used(v_ix)));
    var_used(v_ix) = 1;
    % add the vars to a block
    block = [block, v_ix];    
    b_size = length(block);
    % are we done with the block?
    if(b_size >= block_sizes(cnt_b) || n_group == 1)
      % commit the block
      blocks{cnt_b} = block;
      cnt_b = cnt_b+1;
      block = [];
      n_sum = n_sum + b_size;
    end
    n_group = n_group-1;    
  end

  % sanity check
  if (n_sum ~= n_vars)
    error("make blocks wrong!")
  end

   out.grp_order = [];
   out.blocks = blocks;
end