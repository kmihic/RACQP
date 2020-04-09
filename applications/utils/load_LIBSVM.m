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

function m = load_LIBSVM(filename, p_train, cross_val, fileT, p_test)
rng(1);
  [y_label, x_data] = libsvmread(filename);
  N = length(y_label);
  y_min = min(y_label);
  y_max = max(y_label);
  p_cross_val = 0.3;
  if(y_min ~= -1 || y_max ~= 1)
    y_label(find(y_label == y_min)) = -1;
    y_label(find(y_label == y_max)) = 1;
  end
  if(p_train < 1)
    n_train = ceil(N*p_train);
  else
    n_train = N;
  end
  if(cross_val)
    train_ix = 1:n_train;
    tt_ix = randperm(n_train);
    n_cross = ceil(n_train*p_cross_val);
    test_ix = tt_ix(1:n_cross);
    train_data = x_data(train_ix,:);
    train_label = y_label(train_ix,:);
    test_data = x_data(test_ix,:);
    test_label = y_label(test_ix,:);
  elseif(strcmpi(fileT,"")) %no test file
    tt_ix = randperm(N);
    train_ix = tt_ix(1:n_train);
    test_ix = tt_ix(n_train+1:N);
    train_data = x_data(train_ix,:);
    train_label = y_label(train_ix,:);
    test_data = x_data(test_ix,:);
    test_label = y_label(test_ix,:);
  else
    train_ix = 1:n_train;
    [y_label_t, x_data_t] = libsvmread(fileT); 
    y_min = min(y_label_t);
    y_max = max(y_label_t);
    if(y_min ~= -1 || y_max ~= 1)
     y_label_t(find(y_label_t == y_min)) = -1;
     y_label_t(find(y_label_t == y_max)) = 1;
    end
    N = length(y_label_t);
    if(p_test < 1)
      n_test = ceil(N*p_test);
      tt_ix = randperm(N);
      test_ix = tt_ix(1:n_test);
    else
      test_ix = 1:N;
    end
    train_data = x_data(train_ix,:);
    train_label = y_label(train_ix,:);
    test_data = x_data_t(test_ix,:);
    test_label = y_label_t(test_ix,:);
  end
  
  % construct the output struct
  m = struct();
  m.train_label = train_label;
  m.test_label = test_label;
  m.train_data = train_data;
  m.test_data = test_data;



end






