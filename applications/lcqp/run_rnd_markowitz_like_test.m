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


function run_rnd_markowitz_like_test(solver)
addpath('../racqp');
addpath('../utils');

  inst = get_instances();
  r_time = [];
  n_iter = [];
  p_res = [];
  d_res = [];
  dir = "../data/rnd_markov_like_cont/";
  tmp_dir = "./tmp_rac/";
  mkdir(tmp_dir);
  for jj = 1:length(inst)
    for ii = 1:10
      file_n = inst(jj)+"_RND"+ii;
      filename = file_n+".rac";
      if(length(strfind(filename,"_N6000_")) > 0)
        pos = 1;
      else
        pos = 2;
      end
      disp("Solving: "+file_n)
      disp("Unzipping the file")
      gunzip(dir+filename+".gz", tmp_dir);
      [s, m] = verify_Markowitz_min_variance_rnd(tmp_dir+filename, solver);
      delete(tmp_dir+filename);
      r_time(pos,ii) = s.rac_time;
      n_iter(pos,ii) = s.n_iter;
      p_res(pos,ii) = s.sol_residue_p;
      d_res(pos,ii) = s.sol_residue_d;
    end
  end
  rmdir(tmp_dir,'s')
  disp(" ")
  disp("#####################")
  disp('SUMMARY: Mean values')  
  disp(" ")
  disp("Problem size = 6000");
  disp("Run time: " + mean(r_time(1,:)))
  disp("Num. Iter: " + mean(n_iter(1,:)))
  disp("Primal residual: " + mean(p_res(1,:)))
  disp("Dual residual: " + mean(d_res(1,:)))
  disp(" ")
  disp("Problem size = 9000");
  disp("Run time: " + mean(r_time(2,:)))
  disp("Num. Iter: " + mean(n_iter(2,:)))
  disp("Primal residual: " + mean(p_res(2,:)))
  disp("Dual residual: " + mean(d_res(2,:)))

end

function inst = get_instances()
inst=["MARKOWITZ_N6000_SP0.95_E0.5_Z2_C100";
"MARKOWITZ_N6000_SP0.95_E0.5_Z2_C1.5";
"MARKOWITZ_N9000_SP0.95_E0.5_Z2_C100";
"MARKOWITZ_N9000_SP0.95_E0.5_Z2_C1.5"];

end
