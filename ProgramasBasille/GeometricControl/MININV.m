function Q = mininv(A,B)
% MININV Q = mininv(A,B) is an orthonormal basis for the minimal
%        A-invariant containing imB .
%        The routine implements Algorithm 3.2-1 of "Controlled and
%        Conditioned Invariant Subspaces in Linear System Theory"

%        Basile and Marro 4-20-90

   nv = 0;
   B1 = ima(B,0);
   Q=B1;
   [mq,nq] = size(Q);
   n = 0;
while (nq-nv) > 0
   nv = nq;
   Q = ima([B1 A*Q],0);
   [mq,nq] = size(Q);
% Monitoring iterations
%  n = n+1;
%  disp('   iterations in progress'), n, Q
end
% --- last line of mininv
