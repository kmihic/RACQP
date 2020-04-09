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


function groups = construct_rnd_groups(n_vars, n_blocks)

    if(n_vars < n_blocks)
      error("Number of variables < number of blocks")
    end
    groups={};
    g_size = floor(n_vars/n_blocks);
    n = g_size*n_blocks;
    ix = randperm(n_vars);
    T = reshape(ix(1:n),g_size,[]);
    for ii = 1:size(T,2)
      groups(ii,:) = {T(:,ii)',[]};
    end
    if((n_vars - n) > 0)
      groups(ii+1,:) = {ix(n+1:n_vars),[]};
    end
end
