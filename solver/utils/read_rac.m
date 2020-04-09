
function m = read_rac(filename)

  s = read_rac_file(filename);
  m.size = s.p_size;
  % 0 - x'Qx; 1 - 1/2 x'Qx
  if(s.model_type == 0)
    m.Q = s.Q;
  else
    m.Q = s.Q/2;
  end
  m.c = s.c;
  % constraints
  if(length(s.beq) == 0)
    m.Aeq = zeros(0,m.size);
    m.beq = zeros(0,1);
  else
    m.Aeq = s.Aeq;
    m.beq = s.beq;
  end
  if(length(s.bineq) == 0)
    m.Aineq = zeros(0,m.size);
    m.bineq = zeros(0,1);
  else
    m.Aineq = s.Aineq;
    m.bineq = s.bineq;
  end
  % init x
  m.x0 = zeros(m.size,1);
  % integers
  if(length(s.int_vars) == 0)
    m.integers = [];
    m.binary = [];
  else
    lb_zero = find(~s.lb);
    ub_one = find(s.ub == 1);
    boxed = zeros(m.size,1);
    boxed(lb_zero) = boxed(lb_zero) + 1;
    boxed(ub_one) = boxed(ub_one) + 1;
    boxed_ix = find(boxed == 2);
    all_var = zeros(m.size,1);
    all_var(s.int_vars) = 2; 
    all_var(boxed_ix) = all_var(boxed_ix) + 1;
    % binary are '3', integer '2', cont boxed '1', cont '0'
    m.integers = find(all_var == 2);
    m.binary = find(all_var == 3);
  end

  % bounds
  inf_ix = find(s.lb == -s.inf_val);
  m.lb = s.lb;
  m.lb(inf_ix) = -inf;  
  inf_ix = find(s.ub == s.inf_val);
  m.ub = s.ub;
  m.ub(inf_ix) = inf;
  % constant value
  m.const = 0;
end

function s=read_rac_file(filename, add_to_diagonal)

fileID = fopen(filename,'r');
%    Instance name
s.name = fgetl(fileID);
%    # Using x'Qx or 1/2x'Qx (0/1) model
comment = fgetl(fileID);
s.model_type = sscanf(fgetl(fileID),'%d');
%    # Value for infinity
comment = fgetl(fileID);
s.inf_val = sscanf(fgetl(fileID),'%f');
%    # Problem size
comment = fgetl(fileID);
s.p_size = sscanf(fgetl(fileID),'%d');
%    # Q
comment = fgetl(fileID);
m_size = sscanf(fgetl(fileID),'%d %d %d');
if m_size(3) > 0
  M = fscanf(fileID,'%d %d %f',[3,m_size(3)])';
  Q = get_matrix(M,m_size(1), m_size(1));
  if(nargin>1)
    s.Q = sparse(Q + eye(m_size(1))*add_to_diagonal);
  else
    s.Q = sparse(Q);
  end
  comment = fgetl(fileID); %fscanf does not read newLine in!
else
  s.Q = []; 
end
%    # c
comment = fgetl(fileID);
v_size = sscanf(fgetl(fileID),'%d %d');
if v_size(2) > 0
  V = fscanf(fileID,'%d %f',[2,v_size(2)])';
  s.c = get_vector(V,v_size(1));
  comment = fgetl(fileID);
else
    s.c = [];
end
%    # beq
comment = fgetl(fileID);
v_size = sscanf(fgetl(fileID),'%d %d');
if v_size(2) > 0
  V = fscanf(fileID,'%d %f',[2,v_size(2)])';
  s.beq = get_vector(V,v_size(1));  
  comment = fgetl(fileID);
elseif v_size(1) == 0
        s.beq = [];
else
  s.beq = zeros(v_size(1),1); 
end
%    # Aeq
comment = fgetl(fileID);
m_size = sscanf(fgetl(fileID),'%d %d %d');
if m_size(3) > 0
  M = fscanf(fileID,'%d %d %f',[3,m_size(3)])';
  s.Aeq = sparse(get_matrix(M,m_size(1), m_size(2)));
  comment = fgetl(fileID); %fscanf does not read newLine in!
else
  s.Aeq = []; 
end
%    # bineq
comment = fgetl(fileID);
v_size = sscanf(fgetl(fileID),'%d %d');
if v_size(2) > 0
  V = fscanf(fileID,'%d %f',[2,v_size(2)])';
  s.bineq = get_vector(V,v_size(1));  
  comment = fgetl(fileID);
elseif v_size(1) == 0
        s.bineq = [];
else
  s.bineq = zeros(v_size(1),1); 
end
%    # Aineq
comment = fgetl(fileID);
m_size = sscanf(fgetl(fileID),'%d %d %d');
if m_size(3) > 0
  M = fscanf(fileID,'%d %d %f',[3,m_size(3)])';
  s.Aineq = sparse(get_matrix(M,m_size(1), m_size(2)));
  comment = fgetl(fileID); %fscanf does not read newLine in!
else
  s.Aineq = []; 
end
%    # l_bound
comment = fgetl(fileID);
v_size = sscanf(fgetl(fileID),'%d %d');
if v_size(2) > 0
  V = fscanf(fileID,'%d %f',[2,v_size(2)])';
  s.lb = get_vector(V,v_size(1));  
  comment = fgetl(fileID);
else
  s.lb = zeros(v_size(1),1);  
end
%    # u_bound
comment = fgetl(fileID);
v_size = sscanf(fgetl(fileID),'%d %d');
if v_size(2) > 0
  V = fscanf(fileID,'%d %f',[2,v_size(2)])';
  s.ub = get_vector(V,v_size(1));  
  comment = fgetl(fileID);
else
  s.ub = zeros(v_size(1),1); 
end 
%    # integers
comment = fgetl(fileID);
v_size = sscanf(fgetl(fileID),'%d');
if v_size(1) > 0
  s.int_vars = fscanf(fileID,'%d',[1,v_size(1)])';
  comment = fgetl(fileID);
else
  s.int_vars = []; 
end
fclose(fileID);
end

function K = get_matrix(L,nrows,ncols)
 K = zeros(nrows,ncols);
 for n=1:size(L,1)
  K(L(n,1),L(n,2)) = L(n,3);
 end
end

function K = get_vector(L,nrows)
 K = zeros(nrows,1);
 for n=1:size(L,1)
  K(L(n,1)) = L(n,2);
 end
end


