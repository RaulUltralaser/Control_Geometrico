function Q = ortc(A)
% ORTC	Complementary orthogonalization. Q=ortc(A) is an orthonormal
%	basis for the orthogonal complement of imA .

%	Basile and Marro 4-20-90 (from Matlab's null)

   [Q,R,E] = qr(A);
tol = eps*norm(A,'fro')*10^4;
   [ma,na] = size(A);
if na == 1
   d = R(1,1);
else
   d = diag(R);
end
   d = (abs(d))';
   dd = ma-na;
if dd >= 1
   d=[d zeros(1,dd)];
end
nul = find(d <= tol);
   r=length(nul > 0);
if r > 0
   Q=Q(:,nul);
else
   Q = zeros(ma,1);
end
% --- last line of ortc ---
