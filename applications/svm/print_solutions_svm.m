
function print_solutions_svm(solutions)

  N = length(solutions);
  accuracy = zeros(N,1);
  beta = zeros(N,1);
  sigma = zeros(N,1);
  C = zeros(N,1);
  instance_name = [];
  r_time = [];
  if(isfield(solutions(1),'name'))
    use_name = true;
  else
    use_name = false;
  end
  if(isfield(solutions(1),'beta'))
    use_beta = true;
  else
    use_beta = false;
  end
  for ii = 1:length(solutions)
    accuracy(ii) = solutions(ii).accuracy;
    if(use_beta)
      beta(ii) = solutions(ii).beta;
    end
    C(ii) = solutions(ii).C;
    sigma(ii) = solutions(ii).sigma;
    r_time = [r_time; solutions(ii).r_time];
    if(use_name)
      instance_name = [instance_name; solutions(ii).name];
    end
  end
  T = table(sigma,C,accuracy);
  if(use_beta)
    T = [T,table(beta)];
  end
  if(use_name)
    T = [table(instance_name),T];
  end
  T = [T,struct2table(r_time)];
  disp("###### ")
  disp(" ")
  disp(T);
end
    
