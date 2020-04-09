%  This function is a part of rac_admm_lasso
%  Written for MATLAB_R2016b
%  Copyright (C) 2019
%     Kresimir Mihic <kmihic@alumni.stanford.edu>
%     Mingxi Zhu <mingxiz@stanford.edu>
%     Yinyu Ye <yyye@stanford.edu>

function beta = soft(a,b)
   size = length(b); 
   beta = zeros(size);
   for i = 1:size
    if b(i) > abs(a(i))
        beta(i) = 0;
    else
       if a(i) > 0
           beta(i) = -(a(i)-b(i));
       else
           beta(i) = -(a(i)+b(i));
       end
    end
   end
          
       
end
