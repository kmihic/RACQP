
%
function plot_admm_variants_evol(solutions)

% to make plots, must turn on verbose (run_params.debug=1)
% to produce iter-eval arrays
  rac_p = solutions(3).res_iter;
  rac_d = solutions(3).res_iter_dual;
  rp_p = solutions(6).res_iter;
  rp_d = solutions(6).res_iter_dual;
  cadmm_p = solutions(9).res_iter;
  cadmm_d = solutions(9).res_iter_dual;
  dadmm_p = solutions(12).res_iter;
  dadmm_d = solutions(12).res_iter_dual;

  plot(log(rac_p),'b:','LineWidth',3)
  hold on
  plot(log(rp_p),'r-','LineWidth',3)
  plot(log(cadmm_p),'k-.','LineWidth',3)
  plot(log(dadmm_p),'m--','LineWidth',3)
  axis([0,100,-20,10])
  legend('RAC-ADMM','RP-ADMM','Cyclic-ADMM','Distributed-ADMM')
  xlabel('Iterations')
  ylabel('log(primal residual absolute)')
  title('Primal residual')
  hold off
  figure
  plot(log(rac_d),'b:','LineWidth',3)
  hold on
  plot(log(rp_d),'r-','LineWidth',3)
  plot(log(cadmm_d),'k-.','LineWidth',3)
  plot(log(dadmm_d),'m--','LineWidth',3)
  axis([0,100,-20,10])
  legend('RAC-ADMM','RP-ADMM','Cyclic-ADMM','Distributed-ADMM')
  xlabel('Iterations')
  ylabel('log(dual residual absolute)')
  title('Dual residual')
  hold off
end