
   

function [r_val] = use_libsvm(m, no_test)

  global g_gamma;
  
%   n = length(m.train_label);
%   K_train = [(1:n)', m.K_train];
%   % Precomputed Gaussian Kernel, test
%   n = length(m.test_label);
%   K_test = [(1:n)', m.K_test];
%   n = length(m.train_label);
%   K_train = [(1:n)', m.K_train];
  accuracy = 0;
  
  tstart = tic;
  %opt = sprintf('-t 4 -c %f ',m.C);
  %model = svmtrain(m.train_label, K_train, opt);

  opt = sprintf('-t 2 -c %f -g %f -q ',m.C, g_gamma);
  model_gauss = svmtrain(m.train_label, m.train_data, opt);

  telapsed = toc(tstart);
  disp("LIBSVM test data")
  if(no_test)
    accuracy = -1;
    disp("No testing done");
  else
    [predict_label, accuracy, prob_estimates] = svmpredict(m.test_label, m.test_data, model_gauss);

  end
  disp("Runtime: " + telapsed + " (LIBSVM)")
  r_time = telapsed;
r_val.accuracy = accuracy(1);
r_val.r_time = r_time;
  disp(" ")
end

%accuracy_G

% -s svm_type : set type of SVM (default 0)
%         0 -- C-SVC              (multi-class classification)
%         1 -- nu-SVC             (multi-class classification)
%         2 -- one-class SVM
%         3 -- epsilon-SVR        (regression)
%         4 -- nu-SVR             (regression)
% -t kernel_type : set type of kernel function (default 2)
%         0 -- linear: u'*v
%         1 -- polynomial: (gamma*u'*v + coef0)^degree
%         2 -- radial basis function: exp(-gamma*|u-v|^2)
%         3 -- sigmoid: tanh(gamma*u'*v + coef0)
%         4 -- precomputed kernel (kernel values in training_set_file)
% -d degree : set degree in kernel function (default 3)
% -g gamma : set gamma in kernel function (default 1/num_features)
% -r coef0 : set coef0 in kernel function (default 0)
% -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)
% -n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)
% -p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
% -m cachesize : set cache memory size in MB (default 100)
% -e epsilon : set tolerance of termination criterion (default 0.001)
% -h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)
% -b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)
% -wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)
% -v n: n-fold cross validation mode
% -q : quiet mode (no outputs)




  % Gaussian Kernel
  %g_gamma = 1/(2*sigma^2);
  %opt = sprintf('-t 2 -c %f -g %f ',C, g_gamma);
  %model_gauss = svmtrain(train_label, train_data, opt);
  %[predict_label_G, accuracy_G, prob_estimates] = svmpredict(test_label, test_data, model_gauss);
