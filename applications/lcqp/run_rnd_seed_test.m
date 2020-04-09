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

function run_rnd_seed_test(solver,r_time,epsilon,max_iter)

  rnd_seed = [12 123 1234 12345 123456 654321 54321 4321 321 21];
  solutions = [];
  sol_stat = [];
  inst = get_instances();
  for ii=1:size(inst,1)
    sol = [];
    for rnd_s = rnd_seed
      disp("##### RND SEED: "+rnd_s+"  #####");
      if(strcmpi(inst(ii,1),'markowitz'))
        s = run_portfolio_test(solver,r_time,false,false,epsilon,...
                     max_iter,rnd_s,true,inst(ii,2:end));
      elseif(strcmpi(inst(ii,1),'qap'))
        s = run_qaplib_test(solver,r_time,false,false,epsilon,...
                     max_iter,rnd_s,true,inst(ii,2:end));
      elseif(strcmpi(inst(ii,1),'lptest'))
        s = run_LPTestSet_test(solver,r_time, epsilon, ...
                     max_iter,rnd_s,true,inst(ii,2:end));
      end %if
      s.name = s.name + "_RND" + rnd_s;
      sol = [sol,s];
    end %for rnd seed
    % solution per instance
    solutions = [solutions,sol];
    % stats per instance
    s = get_stats(sol);
    s.name = inst(ii,2);
    sol_stat = [sol_stat,s];
  end %for instance
  disp(" ")
  disp("####################")
  disp(" INFO PER PROBLEM ")
  print_solutions(solutions);    
  disp(" ")
  disp("####################")
  disp(" STATISTICS PER PROBLEM TYPE ")
  disp("####################")    
  print_stats(sol_stat);
end


function inst = get_instances()
  inst=[
  "markowitz", "regular_monthly","","","";
  "markowitz", "regular_daily","","","";
  "qap", "sko100a","","",""; 
  "qap", "tai125e01","","",""; 
  "qap", "tai150b","","",""; 
  "qap", "tho150","","",""; 
  "qap", "wil100","","","";
  "lptest", "square15", "10", "200", "none"];
end





