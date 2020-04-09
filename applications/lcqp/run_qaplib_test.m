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
% Verify results: QAPLIB binary and cont 
%      Using QAPLIB benchmark problems 
%      (http://anjos.mgi.polymtl.ca/qaplib/)
%


function solutions = run_qaplib_test(solver,r_time,binary,all_lib, ...
            epsilon, max_iter, rnd_seed,quiet,inst)

  max_time = r_time;

  if(nargin <= 4)
    epsilon = 1e-6;
  end
  if(nargin <= 5)
    max_iter = 4000;
  end
  if(nargin <= 6)
    rnd_seed = 123;
  end
  if(nargin <= 7)
    quiet = false;
  end
  if(nargin <= 8) 
    if(nargin <= 3 || ~all_lib)
      if(binary)
        inst = get_instances_short();
      else
        inst = get_instances_selected();
      end
    else
        inst = get_instances_long();
    end
  end

  data_path = "../data/data_qaplib/";
  solutions = [];
  for ii = 1:size(inst,1)
    filename = data_path+inst(ii,1);
    disp('Solving: '+filename);
    disp("LOADING THE MODEL...")
    model = load_QAPLIB(filename+'.dat');
    inst_param = get_instance_run_params(model, solver, r_time, binary, ...
                         epsilon, max_iter);
    if(binary)
      inst_param.racqp_run_p = get_qaplib_run_params_bin(model.size, max_time);
    else
      inst_param.racqp_run_p = get_qaplib_run_params_cont(model.size, rnd_seed, ...
                          max_time, epsilon, max_iter, model.Q);
      inst_param.model.binary=[];
    end 
    s = solve_instance(inst_param);
    s.name = inst(ii,1);
    % opt/best obj_val is for binary problems
    if(binary)
      s.obj_val = inst(ii,2); 
    end
    solutions = [solutions,s];
  end
  if(~quiet)
    if(binary) 
      print_solutions_binary(solutions)
    else
      print_solutions(solutions)
    end
  end
end
  
function run_mip = get_qaplib_run_params_bin(nsize,max_time) 

  % MIP run-parameters
  lambda = 0.4;
  p_trial = 0.005;
  %do not split groups
  model.group_mode = 'RAC_NO_SPLIT';
  
  %load default runtime parameters
  r = floor(sqrt(nsize));
  n_blocks = ceil(r/2);
  beta = nsize;
  run_params = default_run_params(n_blocks, beta,max_time);
  %change some parameters
  run_params.max_iter = 1000;
  %turn on verbose
  %run_params.debug = 1;
 
  %load default gurobi parameters
  gurobi_params = default_gurobi_parameters;
  %change how long gurobi will spend on each subproblem
  gurobi_params.TimeLimit = 5;

  run_params.gurobi_params = gurobi_params;

  %set MIP run parameters
  run_mip.permute_type = 'user_defined';
  run_mip.permute_f = @perturb_qap;
  run_mip.permute_dist = 'exponential'; 
  run_mip.permute_min = 2; 
  run_mip.permute_max = r; 
  run_mip.permute_mu = lambda*r; 
  run_mip.max_iter = Inf;
  run_mip.max_nperturb = Inf;
  run_mip.max_rtime = max_time; 
  run_mip.rnd_seed = 123; 
  run_mip.n_perturb_trial = max(2,round(p_trial*nsize));
  run_mip.run_sub = run_params;
  run_mip.debug = 0;
  run_mip.mip_epsilon = 0;
end

function run_params = get_qaplib_run_params_cont(nsize, rnd_seed, max_time, ...
                        epsilon, max_iter, Q) 

  %num blocks = r, block size = r => RP mode
  model.group_mode = 'RP';
  %load default runtime parameters
  r = floor(sqrt(nsize));
  n_blocks = r;
  beta = r;
  run_params = default_run_params(n_blocks, beta,max_time,epsilon,max_iter);
  %change some parameters
  run_params.rnd_seed = rnd_seed;
  %turn on verbose
  %run_params.debug = 1;

  %Check density
  density_Q=nnz(Q)/size(Q,1)^2;
  %model.Aeq is sparse, no need to check
  if(density_Q > 0.5)
    run_params.use_sparse = false;
  else
    run_params.use_sparse = true;
  end
end




function inst = get_instances_selected()
inst=[
%"dre110";
%"esc128";
"dre132";
"sko100a";
%"sko100b";
%"sko100c";
%"sko100d";
%"sko100e";
"sko100f";
"tai100a";
%"tai100b";
"tai125e01";
%"tai150b";
"tho150";
"wil100";
];
end



function inst = get_instances_short()
inst=[
"lipa80a", "253195";
"lipa80b", "7763962";
"lipa90a", "360630";
"lipa90b", "12490441";
"sko81", "90998";
"sko90", "115534";
"sko100a", "152002";
"sko100b", "153890";
"sko100c", "147862";
"sko100d", "149576";
"sko100e", "149150";
"sko100f", "149036";
"tai80a", "13499184";
"tai80b", "818415043";
"tai100a", "21052466";
"tai100b", "1185996137";
"tai150b", "498896643";
"tho40", "240516";
"tho150", "8133398";
"wil50", "48816";
"wil100", "273038";
];
end

function inst = get_instances_long()
inst=[
"bur26a", "5426670";
"bur26b", "3817852";
"bur26c", "5426795";
"bur26d", "3821225";
"bur26e", "5386879";
"bur26f", "3782044";
"bur26g", "10117172";
"bur26h", "7098658";
"chr12a", "9552";
"chr12b", "9742";
"chr12c", "11156";
"chr15a", "9896";
"chr15b", "7990";
"chr15c", "9504";
"chr18a", "11098";
"chr18b", "1534";
"chr20a", "2192";
"chr20b", "2298";
"chr20c", "14142";
"chr22a", "6156";
"chr22b", "6194";
"chr25a", "3796";
"els19", "17212548";
"esc128", "64";
"esc16a", "68";
"esc16b", "292";
"esc16c", "160";
"esc16d", "16";
"esc16e", "28";
"esc16g", "26";
"esc16h", "996";
"esc16i", "14";
"esc16j", "8";
"esc32a", "130";
"esc32b", "168";
"esc32c", "642";
"esc32d", "200";
"esc32e", "2";
"esc32g", "6";
"esc32h", "438";
"esc64a", "116";
"had12", "1652";
"had14", "2724";
"had16", "3720";
"had18", "5358";
"had20", "6922";
"kra30a", "88900";
"kra30b", "91420";
"kra32", "88700";
"lipa20a", "3683";
"lipa20b", "27076";
"lipa30a", "13178";
"lipa30b", "151426";
"lipa40a", "31538";
"lipa40b", "476581";
"lipa50a", "62093";
"lipa50b", "1210244";
"lipa60a", "107218";
"lipa60b", "2520135";
"lipa70a", "169755";
"lipa70b", "4603200";
"lipa80a", "253195";
"lipa80b", "7763962";
"lipa90a", "360630";
"lipa90b", "12490441";
"nug12", "578";
"nug14", "1014";
"nug15", "1150";
"nug16a", "1610";
"nug16b", "1240";
"nug17", "1732";
"nug18", "1930";
"nug20", "2570";
"nug21", "2438";
"nug22", "3596";
"nug24", "3488";
"nug25", "3744";
"nug27", "5234";
"nug28", "5166";
"nug30", "6124";
"rou12", "235528";
"rou15", "354210";
"rou20", "725522";
"scr12", "31410";
"scr15", "51140";
"scr20", "110030";
"sko100a", "152002";
"sko100b", "153890";
"sko100c", "147862";
"sko100d", "149576";
"sko100e", "149150";
"sko100f", "149036";
"sko42", "15812";
"sko49", "23386";
"sko56", "34458";
"sko64", "48498";
"sko72", "66256";
"sko81", "90998";
"sko90", "115534";
"ste36a", "9526";
"ste36b", "15852";
"ste36c", "8239110";
"tai100a", "21043560";
"tai100b", "1185996137";
"tai12a", "224416";
"tai12b", "39464925";
"tai150b", "498896643";
"tai15a", "388214";
"tai15b", "51765268";
"tai17a", "491812";
"tai20a", "703482";
"tai20b", "122455319";
"tai256c", "44759294";
"tai25a", "1167256";
"tai25b", "344355646";
"tai30a", "1818146";
"tai30b", "637117113";
"tai35a", "2422002";
"tai35b", "283315445";
"tai40a", "3139370";
"tai40b", "637250948";
"tai50a", "4938796";
"tai50b", "458821517";
"tai60a", "7205962";
"tai60b", "608215054";
"tai64c", "1855928";
"tai80a", "13499184";
"tai80b", "818415043";
"tho150", "8133398";
"tho30", "149936";
"tho40", "240516";
"wil100", "273038";
"wil50", "48816";
];
end
