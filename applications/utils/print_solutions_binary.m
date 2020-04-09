

function print_solutions_binary(solutions, first_col_title,msg)

  if(nargin <= 1)
   first_col_title = 'Instance_name';
  end
  if(nargin <= 2)
    msg = "";
  end

  name = [];
  note = [];
  N = length(solutions);
  obj_val = zeros(N,1);
  best_obj_val = zeros(N,1);
  gap = zeros(N,1);
  run_time = zeros(N,1);
  res_p = zeros(N,1);
  for ii = 1: length(solutions)
    name = [name;solutions(ii).name];
    obj_val(ii) = solutions(ii).sol_obj_val;    
    best_obj_val(ii) = solutions(ii).obj_val;
    gap(ii) = (obj_val(ii)-best_obj_val(ii))/(1+abs(best_obj_val(ii)));
    run_time(ii) = solutions(ii).runtime;
    res_p(ii) = solutions(ii).sol_res_p.rel;
    if(~isfinite(solutions(ii).sol_res_p.rel))
      note = [note;string("Time expired. Heuristic/best found solution")];
    elseif(solutions(ii).sol_res_p.rel ~= 0)
      note = [note;string("Primal residual: "+ solutions(ii).sol_res_p.rel)];
    else
      note = [note; string("--")];
    end
  end
  T = table(name, run_time, best_obj_val,obj_val, gap,note);  
  T.Properties.VariableNames={first_col_title,'run_time','opt_best_known_obj_val',...
      'obj_val','gap','note'};  
  disp(" ")
  disp("####")
  disp(msg);
  disp(" ");
  disp(T);
end