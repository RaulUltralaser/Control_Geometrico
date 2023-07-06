function R = mlocus(A,B,C,D,k)
%RLOCUS Multivariable root locus.
%       R = MLOCUS(A,B,C,D,k) finds the root locus for the multivariable
%       system:
%               .
%               x = Ax + Bu    u = -k*y
%               y = Cx + Du
%
%       The number of outputs must be equal to that of inputs.
%       MLOCUS returns a matrix R with LENGTH(k) rows and as many columns
%       as the state space dimension is (dimension of the square matrix A).
%       Each row of the matrix corresponds to a gain from vector k.
%       The root locus may be plotted with  PLOT(R,'x').
%       The program is a generalized version of Matlab's rlocus.m

%       Basile and Marro 4-20-90

[t,na] = size(A);
[t,nb] = size(B);
[mc,t] = size(C);
if nb ~= mc
  error('   different number of inputs and outputs in MLOCUS');
end
error(abcdchk(A,B,C,D));
nk = length(k);
R  = sqrt(-1) * ones(na,nk);
% Find eigenvalues of:  A - B*inv(I+k*D)*k*C:
Id = eye(nb);
for i=1:nk
  R(:,i) = eig(A-k(i)*B*inv(Id+k(i)*D)*C);
end
R = R.';
% --- last line of mlocus ---
