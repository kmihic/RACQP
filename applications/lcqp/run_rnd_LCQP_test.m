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


function run_rnd_LCQP_test(solver)
addpath('../racqp');
addpath('../utils');

  inst = get_instances();
  r_time = [];
  n_iter = [];
  p_res = [];
  d_res = [];
  dir = "../data/rnd_lcqp/";
  tmp_dir = "./tmp_rac/";
  mkdir(tmp_dir,'s');
  for jj = 1:length(inst)
    for ii = 1:10
      file_n = inst(jj)+"_RND"+ii;
      filename = file_n+".rac";
      if(length(strfind(filename,"_N6000_")) > 0)
        if(length(strfind(filename,"_M600_")) > 0)
          pos = 1;
        else
          pos = 2;
        end
      else
        if(length(strfind(filename,"_M900_")) > 0)
          pos = 3;
        else
          pos = 4;
        end
      end
      disp("Solving: "+file_n)
      disp("Unzipping the file")
      gunzip(dir+filename+".gz", tmp_dir);
      [s, m] = verify_rnd_LCQP(tmp_dir+filename, solver);
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
  disp("Number of rows = 600 (for each Aeq, Aineq)")
  disp("Run time: " + mean(r_time(1,:)))
  disp("Num. Iter: " + mean(n_iter(1,:)))
  disp("Primal residual: " + mean(p_res(1,:)))
  disp("Dual residual: " + mean(d_res(1,:)))
  disp(" ")
  disp("Number of rows = 3000 (for each Aeq, Aineq)")
  disp("Run time: " + mean(r_time(2,:)))
  disp("Num. Iter: " + mean(n_iter(2,:)))
  disp("Primal residual: " + mean(p_res(2,:)))
  disp("Dual residual: " + mean(d_res(2,:)))
  disp(" ")
  disp("Problem size = 9000");
  disp("Number of rows = 900 (for each Aeq, Aineq)")
  disp("Run time: " + mean(r_time(3,:)))
  disp("Num. Iter: " + mean(n_iter(3,:)))
  disp("Primal residual: " + mean(p_res(3,:)))
  disp("Dual residual: " + mean(d_res(3,:)))
  disp(" ")
  disp("Problem size = 9000");
  disp("Number of rows = 4500 (for each Aeq, Aineq)")
  disp("Run time: " + mean(r_time(4,:)))
  disp("Num. Iter: " + mean(n_iter(4,:)))
  disp("Primal residual: " + mean(p_res(4,:)))
  disp("Dual residual: " + mean(d_res(4,:)))
end

function inst = get_instances()
inst=["LCQP_N6000_SP0.95_M600_E0.5_Z2_C100";
"LCQP_N6000_SP0.95_M600_E0.5_Z2_C1.5";
"LCQP_N6000_SP0.95_M3000_E0.5_Z2_C100";
"LCQP_N6000_SP0.95_M3000_E0.5_Z2_C1.5";
"LCQP_N9000_SP0.95_M900_E0.5_Z2_C100";
"LCQP_N9000_SP0.95_M900_E0.5_Z2_C1.5";
"LCQP_N9000_SP0.95_M4500_E0.5_Z2_C100";
"LCQP_N9000_SP0.95_M4500_E0.5_Z2_C1.5"];

end
