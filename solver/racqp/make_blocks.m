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

% Must make sure that all the vars are specified by the groups struct
% i.e. a single var is also a group by itself, e.g {[3] []}
% Note that by using overlap field, we can construct less random way
% of selecting vars. Consider an example of linking groups in a circular way:
% if overlap field only specifies one group, we end up with same groups,
% following the same sequence of optimization with the only difference in 
% first/last group.

% rp = true -> one group, one block
% no_split = true -> a group must not be split between blocks 
function out = make_blocks(n_vars,groups, block_sizes, no_split)

  n_group = size(groups,1);

  b_ix = 1;
  block = [];  
  grp_order = [];
  grp_cnt=1;
  
  rp = (length(block_sizes) == 0);

  var_used = zeros(n_vars,1);
  grp_used = zeros(n_group);
  grp = randi(n_group); 
  if(rp)
    b_size = length(groups{grp,1});
    n_blocks = inf;
    blocks = cell(1,1);
  else
    b_size = block_sizes(b_ix);
    n_blocks = length(block_sizes);
    blocks = cell(1,1);
  end
  b_map = [1:n_group];
  n_sum = 0;
  while (true) 
    while(b_ix <= n_blocks && length(groups)>0)
      space_left = b_size - length(block);
      v_ix = find(var_used(groups{grp,1})==0);
      v_r = find(var_used(groups{grp,1})==1);
      grp_size = length(v_ix); 
      if(grp_size == 0) %all grp vars already added somewhere
        groups{grp,1}(:)=[]; %remove the group vars
        break;
      end
      if(space_left >= grp_size || no_split)
        block = [block, groups{grp,1}(v_ix)]; %add all vars
        var_used(groups{grp,1}(v_ix)) = 1;
        groups{grp,1}(:)=[];   %remove the group vars          
        n_sum = n_sum + length(v_ix);
        %grp is empty, but block may not be full
      else      
        groups{grp,1}(v_r)=[]; %remove used   
        if(grp_size == 1)
          x_ix = 1;
        else
          x_ix = datasample([1:grp_size],space_left,'Replace',false);
        end
        block = [block, groups{grp,1}(x_ix)]; %add chosen vars
        var_used(groups{grp,1}(x_ix)) = 1;           
        n_sum = n_sum + length(x_ix);
        %block is full, but the group is not empty
      end
      if(length(block) >= b_size || n_sum == n_vars)
          blocks{b_ix} = block;
          if(n_sum == n_vars)  %we are done
            break;
          end
          b_ix = b_ix + 1;
          block = [];
          if(~rp)
            b_size = block_sizes(b_ix);
          end
      end    
    end    
    if(length(b_map) == 0)
      break;
    end
    grp_order(grp_cnt)=b_map(grp);
    %we are done
    if(n_sum == n_vars)  
      break;
    end
    %remove the group from the cell array      
    grp_cnt = grp_cnt+1;
    o_blocks = groups{grp,2};
    groups(grp,:)=[];
    b_map(grp)=[];
    s = size(groups,1); 
    if(s>0) 
      grp = randi(s);      
      if(rp)
        b_size = length(groups{grp,1});
      end
    else
     % no more groups, save the last one
     if(length(block) > 0)
       blocks{b_ix} = block;
     end
    end

  end

  out.blocks = blocks;
  out.grp_order = grp_order;
  
%   if(n_sum ~= n_vars)    
%     error('Error. Grouping not done correctly. Not all variables used.')
%   end

end