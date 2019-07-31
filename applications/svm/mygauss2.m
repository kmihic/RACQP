function K_rbf = mygauss2(u)

  global g_gamma;
  
  utu = u*u';
  du = diag(utu);
  
  K_rbf=exp(-g_gamma*(du-2*utu+du'));
end

