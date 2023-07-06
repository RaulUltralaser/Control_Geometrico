function z = gazero(A,B,C,D)
% GAZERO Invariant zeros.
%        z = GAZERO(A,B,C) or z = GAZERO(A,B,C,D) returns in a column vector
%        the invariant zeros of the state-space system
%               .
%               x = Ax + Bu
%               y = Cx (+ Du)
%
%        The algorithm is based on the standard geometric approach definition
%        of invariant zeros. A different algorithm is used in the program
%        tzero of the standard Matlab Control System Toolbox .

%        Basile and Marro 4-20-90

tol = norm(A,'fro')*eps*100;
error(nargchk(3,4,nargin));
[ma,na] = size(A);
[mb,nb] = size(B);
[mc,nc] = size(C);
k1 = 1;
nargs = nargin;
if nargs == 3
  k1 = 0;
  D=zeros(mc,nb);
end
if norm(D,'fro') < tol
  k1 = 0;
end
error(abcdchk(A,B,C,D));
if(k1 == 1)
  A=[A zeros(na,mc); C zeros(mc)];
  B=[B; D];
  C=[zeros(mc,na) eye(mc)];
end
E = ORTC(C');
J = MAINCO(A,B,E);
[P,Q] = STABV(A,B,J);
z = eig(P);
% --- last line of gazero ---
