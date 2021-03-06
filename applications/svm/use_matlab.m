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


function [r_val] = use_matlab(m, no_test)

  tstart = tic;
  model = fitcsvm(full(m.train_data),full(m.train_label),'KernelFunction','mygauss_mat','BoxConstraint',m.C,'ClassNames',[-1,1]);
  % Matlab uses linear for two-class classification. Cannot force it to do gauss!?!
  telapsed = toc(tstart);
  if(no_test)
    accuracy = -1;
  else
    % get prediction
    y_alpha = zeros(model.NumObservations,1);
    y_alpha(model.IsSupportVector) = model.Alpha;
    y_alpha = y_alpha .* m.train_label;

    label = classify(m.test_data, m.train_data, y_alpha, model.Bias);

    correct = find(~(label - m.test_label));
    nc = length(correct);
    na = length(m.test_label);
    accuracy = 100*nc/na;
  
    disp("Matlab test data")
    msg = sprintf('Accuracy = %.4f%% (%d/%d) (MAT)',accuracy,nc, na);
    disp(msg)
  end
  disp("Runtime: " + telapsed + " (MAT)")
  r_time = telapsed;
  r_val.accuracy = accuracy;  
  r_val.r_time.runtime = r_time;
  disp(" ")
end

