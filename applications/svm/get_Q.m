function Q = get_Q(i,j,x1,x2)

 global Y_train;
 global X_train;
 global X_train_sq;
 %global g_gamma;

  % x'Qx used only for objVal calc. Not needed for SVM
  if(nargin == 4)
    Q = NaN;
    return;
  end
 
 K = bsxfun(@minus,X_train_sq(i),(2*X_train(i,:))*X_train(j,:).');
 K = bsxfun(@plus,X_train_sq(j).',K);

 Q = exp(-K);
 if(nargin == 3)
   Q = Q*x1;
 end

end

 
