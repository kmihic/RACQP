function [free avaiable] = get_free_memory(use_swap)

  [r,w] = unix('free -b | grep Mem');
  stats = str2double(regexp(w, '[0-9]*', 'match'));
  free = stats(3);
  if(length(stats) == 6)
    avaiable = stats(6);
  else
    % this is not so good estimate as it does not take
    % into account page cache and also that not all 
    % reclaimable memory slabs will be reclaimed due to items being in use
    % see "man free"
    avaiable = stats(3)+stats(end);
  end
  if(nargin == 1 && use_swap == true)
    [r,w] = unix('free -b | grep Swap');
    stats = str2double(regexp(w, '[0-9]*', 'match'));
    avaiable = avaiable + stats(3);
    free = free + stats(3);
  end
end