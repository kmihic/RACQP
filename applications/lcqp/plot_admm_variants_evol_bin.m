
%
function plot_admm_variants_evol_bin(solutions)

  rac_p = solutions(1).obj_iter;
  rp_p = solutions(2).obj_iter;
  cadmm_p = solutions(3).obj_iter;
  %dadmm_p = solutions(4).obj_iter;

  figure
  plot(rac_p,'b:','LineWidth',3)
  hold on
  plot(rp_p,'r-','LineWidth',3)
  plot(cadmm_p,'k-.','LineWidth',3)
 % plot(dadmm_p,'m--','LineWidth',3)
  legend('RAC-ADMM','RP-ADMM','Cyclic-ADMM');%,'Distributed-ADMM')
  xlabel('Iterations')
  ylabel('objective value')
  hold off
end
