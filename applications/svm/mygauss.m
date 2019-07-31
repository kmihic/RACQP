function K_rbf = mygauss(u,v)

  global g_gamma;

  u_sq=sum(u.^2,2);
  v_sq=sum(v.^2,2);
  K=bsxfun(@minus,u_sq,(2*u)*v.');
  K=bsxfun(@plus,v_sq.',K);
  
  K_rbf=exp((-g_gamma)*K);
end
