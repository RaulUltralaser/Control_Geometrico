function V = compl(A,B,E,Pl,in)
% COMPL Complementation. Let imVm be the maximal (A,B)-controlled
%       invariant contained in kerE such that sums(imVm,imPl) = X .
%       V = compl(A,B,E,Pl,in) provides an orthonormal V such that
%       sums(imV,imPl) = X, ints(imV,imPl) = {0}, imV contained in kerE .
%       If constant in is equal to 1 some information on the performed
%       change of basis is displayed .
%       The routine implements Algorithm 6.2-1 of "Controlled and
%       Conditioned Invariant Subspaces in Linear System Theory" .

%       Basile and Marro 4-20-90

xx = 'strike any key to continue';
[t,na] = size(A);
tol = norm(A,'fro')*eps*10^4;
kerE = ortc(E');
Vm = mainco(A,B,kerE);
[t,ne] = size(sums(Vm,Pl));
if ne < na
  error('   complementation impossible in COMPL');
end
T = ints(Vm,Pl);
[t,ne] = size(T);
no = norm(T,'fro');
if (ne == 1)&(no < tol)
  disp('   complementation not necessary in compl');
  V = Vm;
  return
end
mic = miinco(A,kerE,B);
Rv = ints(Vm,mic);
T = ima(Rv,0);
[t,n1] = size(T);
no = norm(T,'fro');
if (n1 == 1)&(no < tol)
  n1 = 0;
end
Rp = ints(Vm,Pl);
if n1 == 0
  T = ima(Rp,0);
  Dmic = mic;
  Dmicv = Dmic;
  [t,nz1] = size(Dmic);
else
  Tx = T;
  T = ima([Tx Rp],0);
  Dmicv = ima([Tx mic],0);
  [t,nz1] = size(Dmicv);
  Dmic = Dmicv(:,(n1+1):nz1);
end
  [t,n2] = size(T);
  no = norm(T,'fro');
  if (n2 == 1)&(no < tol)
  n2 = 0;
end
if n2 == 0
  T = ima(Vm,0);
  Dpla = ima(Pl,0);
else
  Tx = T;
  T = ima([Tx Vm],0);
  Dpla = ima([Tx Dmic Pl],0);
  [t,nz2] = size(Dpla);
  Dpla = Dpla(:,(n2+nz1+1):nz2);
end
[t,n3] = size(T);
T = ima([T eye(na)],0);
[t,nz1] = size(Dmic);
T(:,(n3+1):(n3+nz1)) = Dmic(:,1:nz1);
n4 = n3+nz1;
[t,nz2] = size(Dpla);
T(:,(n4+1):(n4+nz2)) = Dpla(:,1:nz2);
if in == 1
  disp('   transformation matrix in compl')
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
  T3 = T(:,1:n3);
  T3
  disp(xx)
  pause
  disp(' ')
  T
  disp(xx)
  pause
  disp(' ')
end
At = inv(T)*A*T;
if in == 1
  disp('   transformed matrices')
  At
  disp(xx)
  pause
  disp(' ')
  Bt = inv(T)*B;
  Bt
  disp(xx)
  pause
  disp(' ')
  Et = E*T;
  Et
  disp(xx)
  pause
  disp(' ')
end
if n2 == n1
  V=ima(T(:,(n2+1):n3),1);
else
  A22 = At((n1+1):n2,(n1+1):n2);
  A33 = At((n2+1):n3,(n2+1):n3);
  A23 = At((n1+1):n2,(n2+1):n3);
  X = lyap(A22,-A33,A23);
  V = ima(T(:,(n1+1):n2)*X+T(:,(n2+1):n3),1);
  if in == 1
    disp('   submatrices for the Sylvester equation')
    A22
    disp(xx)
    pause
    disp(' ')
    A33
    disp(xx)
    pause
    disp(' ')
    A23
    disp(xx)
    pause
    disp(' ')
    disp('   solution of the Sylvester equation')
    X
    disp(xx)
    pause
    disp(' ')
  end
end
%[t,nv] = size(V);
%V = V*inv(V((na-nv+1):na,:));
% --- last line of compl ---
