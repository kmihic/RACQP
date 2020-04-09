
function [A b]= remove_empty_rows(A, b)
  
  B = A;
  B(B~=0)=1;
  s = sum(B,2);
  s_ix = find(s);
  A = A(s_ix,:);
  b = b(s_ix);
end