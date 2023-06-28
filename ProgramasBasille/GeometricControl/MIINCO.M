function Q = miinco(A,C,X)
% MIINCO Q = miinco(A,C,X) is an orthonormal basis for the minimal
%        (A,imC)-conditioned invariant containing imX .
%        The routine implements Algorithm 4.1-1 of "Controlled and
%        Conditioned Invariant Subspaces in Linear System Theory"

%        Basile and Marro 4-20-90

   nv = 0;
   X1 = ima(X,0);
   Q=X1;
   [mq,nq] = size(Q);
   n = 0;
while (nq-nv) > 0
   nv = nq;
   Q = ima([X1 A*(ints(Q,C))],0);
   [mq,nq] = size(Q);
% Monitoring iterations
%  n = n+1;
%  disp('   iterations in progress'), n, Q
end
% --- last line of miinco ---
