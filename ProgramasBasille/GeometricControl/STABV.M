function [P,Q] = stabv(A,B,X)
% STABV  [P,Q] = stabv(A,B,X) provides as P and Q the critical matrices
%	 for the internal and external stabilizability of an (A,B)-
%	 controlled invariant X .

%	 Basile and Marro 4-20-90

% Check and message
   [mx,nx] = size(X);
   T=mainco(A,B,X);
   [my,ny] = size(T);
if (ny ~= nx)
   error('   X is not a controlled invariant in STABV')
end
%
   X = T;
   tol = norm(A,'fro')*eps*10^4;
   [ma,na] = size(A);
   P = miinco(A,X,B);
   R = ints(X,P);
   [m1,n1] = size(R);
   no = norm(R,'fro');
if (n1 == 1)&(no < tol)
   n1 = 0;
end
if n1 == 0
   T = ima(X,0);
   Q = P;
else
   T = ima([R X],0);
   Q = ima([R P],0);
end
   [m2,n2] = size(T);
   [mx,nx] = size(Q);
   Q = Q(:,(n1+1):nx);
   no = norm(T,'fro');
if (n2 == 1)&(no < tol)
   n2 = 0;
end
if n2 == 0
   T = ima([P ortc(P)],0);
else
   T = ima([T ortc(T)],0);
   T(:,(n2+1):(n2+nx-n1)) = Q;
end
   [m3,n3] = size(T);
if (n3 == 1)&(no < tol)
   n3 = 0;
   T = eye(m3);
end
   A1 = inv(T)*A*T;
   ni = n1+1;
   nx = n2-n1;
   P = A1(ni:nx,ni:nx);
   [T Q] = stabi(A,sums(X,mininv(A,B)));
% --- last line of stabv ---
