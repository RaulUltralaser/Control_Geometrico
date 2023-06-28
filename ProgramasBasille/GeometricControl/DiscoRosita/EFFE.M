function F = effe(A,B,X)
% EFFE   F = effe(A,B,X) is a state-to-input feedback matrix such that
%        the (A,imB)-controlled invariant X is an (A+B*F)-invariant .
%        The routine implements Algorithm 4.1-3 of "Controlled and
%        Conditioned Invariant Subspaces in Linear System Theory"

%        Basile and Marro 4-20-90

% Checks and messages
   [mx,nx] = size(X);
   [my,ny] = size(mainco(A,B,X));
if ny ~= nx
error('   X is not a controlled invariant in EFFE')
end
   [mb,nb] = size(B);
   [my,ny] = size(ima(B,1));
if ny ~= nb
error('   B is not full rank in EFFE')
end
%
   tol = norm(A,'fro')*eps*10^4;
   [ma,na] = size(A);
   T = ints(X,B);
   T = ima(T,1);
   [m1,n1] = size(T);
   no = norm(T,'fro');
if (n1 == 1)&(no < tol)
   n1 = 0;
end
if n1 == 0
   T = ima(X,1);
else
   T = ima([T X],0);
end
   [m2,n2] = size(T);
   no = norm(T,'fro');
if (n2 == 1)&(no < tol)
   n2 = 0;
end
if n2 == 0
   T = ima(B,1);
else
   T = ima([T B],0);
end
   [m3,n3] = size(T);
   no = norm(T,'fro');
if (n3 == 1)&(no < tol)
   n3 =0;
   T = eye(m3);
else
   T = ima([T ortc(T)],0);
end
   c1 = n1+1;
   c2 = n3+1;
   T = [B T(:,c1:n2) T(:,c2:na)];
   A1 = inv(T)*A*T;
   c1 = na - n3;
   F = [-A1(1:nb,1:n3) zeros(nb,c1)]*inv(T);
% --- last line of effe ---
