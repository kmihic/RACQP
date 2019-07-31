function m = prepare_data(filename, p_train, fileT, fileT_p)

  [y_label, x_data] = libsvmread(filename);
  N = length(y_label);
  y_min = min(y_label);
  y_max = max(y_label);
 if(y_min ~= -1 || y_max ~= 1)
    y_label(find(y_label == y_min)) = -1;
    y_label(find(y_label == y_max)) = 1;
  end
  n_train = ceil(N*p_train);
  % Split Data
  train_data = x_data(1:n_train,:);
  train_label = y_label(1:n_train,:);
  if(strcmpi(fileT,""))
    be = n_train+1;
    test_size = N-n_train;
    test_size = ceil(test_size*fileT_p);
    bs = min(N,be+test_size);
    test_data = x_data(be:bs,:);
    test_label = y_label(be:bs,:);
  else
    [y_label, x_data] = libsvmread(fileT); 
    y_min = min(y_label);
    y_max = max(y_label);
    if(y_min ~= -1 || y_max ~= 1)
     y_label(find(y_label == y_min)) = -1;
     y_label(find(y_label == y_max)) = 1;
    end
    n_test = ceil(length(y_label)*fileT_p);
    test_data = x_data(1:n_test,:);
    test_label = y_label(1:n_test,:);
  end
  % dataset issue with LIBSVM benchamrs, train features != test ?!?
  n_f_train = size(train_data,2);
  n_f_test = size(test_data,2);
  if(n_f_train < n_f_test)
    test_data(:,n_f_train+1:n_f_test)=[];
  elseif(n_f_train > n_f_test)
    train_data(:,n_f_test+1:n_f_train)=[];
  end

  
  % construct the output struct
  m.train_label = train_label;
  m.test_label = test_label;
  m.train_data = train_data;
  m.test_data = test_data;



end






