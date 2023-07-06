function Q = ima(A,p)
% IMA	Orthogonalization. Q=ima(A,p) is an orthonormal basis for imA .
%	If p=1 permutations are allowed, while if p=0 they are not .

%	Basile and Marro 4-20-90

tol = eps*norm(A,'fro')*10^4;
[ma,na] = size(A);
if p == 1
[Q,R,E] = qr(A);
  if (na == 1)|(ma == 1)
    d = R(1,1);
  else
    d = diag(R);
  end
  d = (abs(d))';
  nul = find(d > tol);
  r=length(nul > 0);
  if r > 0
    Q = Q(:,nul);
  else
    Q = zeros(ma,1);
  end
else
  ki = 1;
  A1=A;
  while ki == 1
    [ma,na] = size(A1);
    punt = 1:na;
    [Q,R] = qr(A1);
    if (na == 1)|(ma == 1)
      d = R(1,1);
    else
      d = diag(R);
    end
    d = (abs(d))';
    nul = find(d <= tol);
    n = min(nul);
    punt = find(punt ~= n);
    if (length(nul) == 0)|(length(punt) == 0)
      ki = 0;
    else
      A1 = A1(:,punt);
    end
  end
  if (length(n) == 1)&(length(punt) == 0)
    Q = zeros(ma,1);
  else
    r = length(d);
    if r > 0
      Q=Q(:,1:r);
      Q = -Q;
      Q(:,r) = -Q(:,r);
    end
  end
end
% --- last line of ima ---
