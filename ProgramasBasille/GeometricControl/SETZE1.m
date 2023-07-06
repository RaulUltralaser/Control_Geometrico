function [B,C,D] = setze1(A,B,C,D,P)
% SETZE  Set zeros.
%	 [B,C,D] = SETZE1(A,B,C,D) or [B,C,D] = SETZE1(A,B,C,D,P) . If called
%	 with D = O it sets C and D (if nb < mc) or B and D (if nb >= mc)
%	 such that the transmission zeros of the system
%		.
%		x = Ax + Bu
%		y = Cx + Du
%
%	 are those communicated in interactive mode or specified in P .
%	 If called with D ~= O, it maintains the communicated value for D .
%	 The routine implements Algorithm 4.5-1 of "Controlled and
%	 Conditioned Invariant Subspaces in Linear System Theory"

%	 Basile and Marro 8-10-90

nargs = nargin;
if (nargs ~= 4)&(nargs ~= 5)
  error('   wrong number of arguments in SETZE1')
end
[t,na] = size(A);
[t,nb] = size(B);
[mc,t] = size(C);
s = 1;
if nb < mc
  s = 0;
end
if s == 0
  A1 = A;
  B1 = B;
  C1 = C;
  D1 = D;
else
  A1 = A';
  B1 = C';
  C1 = B';
  D1 = D';
  r = nb;
  nb = mc;
  mc = r;
end
if nargs == 4
  F1 = place1(A1,B1);
else
  F1 = place1(A1,B1,P);
end
if norm(D,'fro') == 0
  D1 = zeros(mc,nb);
  j = 1;
  for i = 1:mc
    D1(i,j) = 1;
    if j == nb
      j = j - nb;
    end
    j = j+1;
  end
end
if s == 0
  C = -D1*F1;
  D = D1;
else
  B = -F1'*D1';
  D = D1';
end
% --- last line of setze1 ---
