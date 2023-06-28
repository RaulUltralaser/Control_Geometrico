function F = stabf(A,B,V,Pl,in,P)
% STABF State feedback matrix . Let imV be an (A,B)-controlled invariant
%       satisfying sums(imV,imPl) = X, ints(imV,imPl) = 0, and (A1,B1) be
%       stabilizable (A1 is formed with the first n1 rows and columns of A,
%       B1 with the first n1 rows of B: n1 is the number of columns of Pl) .
%       F = staf(A,B,V,Pl,in) or F = staf(A,B,V,Pl,P,in) returns a state
%       feedback matrix such that (A+BF)imV is contained in imV and A1+B1F1
%       has all free eigenvalues assigned . When called with five arguments,
%       the eigenvalues to be assigned are communicated in interactive mode,
%       when called with six arguments they are communicated in column
%       vector P . If constant in is equal to 1 some information on the
%       performed change of basis is displayed .
%       The routine implements Algorithm 6.2-2 of "Controlled and
%       Conditioned Invariant Subspaces in Linear System Theory" .

%       Basile and Marro 4-20-90

[t,na] = size(A);
[t,nb] = size(B);
% Checks and messages
  [t,nv] = size(V);
  [t,ny] = size(mainco(A,B,V));
if ny ~= nv
  error('   V is not a controlled invariant in STABF');
end
  [t,ny] = size(ima(B,1));
if ny ~= nb
  error('   B is not full rank in STABF');
end
xx = 'strike any key to continue';
T = ima(B,1);
[t,n1] = size(T);
T = ima([T Pl],0);
[t,n2] = size(T);
T = ima([T V],0);
[t,n3] = size(T);
if n3 ~= na
  error('   loss of rank in STABF');
end
T(:,1:n1) = B;
T(:,(n2+1):n3) = V;
if in == 1
  disp(' ')
  disp('   transformation matrix in stabf')
  n1
  n2
  n3
  disp(xx)
  pause
  disp(' ')
  T1 = T(:,1:n1);
  T1
  disp(xx)
  pause
  disp(' ')
  T2 = T(:,1:n2);
  T2
  disp(xx)
  pause
  disp(' ')
  T
  disp(xx)
  pause
  disp(' ')
end
At = inv(T)*A*T;
Bt = inv(T)*B;
nargs = nargin;
if nargs == 5
  j=sqrt(-1);
  F1t = place1(At(1:n2,1:n2),Bt(1:n2,:));
else
  [t,np] = size(P);
  if np == 0
    F1t = zeros(nb,n2);
  else
    F1t = place1(At(1:n2,1:n2),Bt(1:n2,:),P);
  end
end
F2t = -At(1:n1,(n2+1):na);
Ft = [F1t F2t];
if in == 1
  disp('   transformed matrices')
  At
  disp(xx)
  pause
  disp(' ')
  Bt
  disp(xx)
  pause
  disp(' ')
  Ft
  disp(xx)
  pause
  disp(' ')
end
F = Ft*inv(T);
% --- last line of stabf ---
