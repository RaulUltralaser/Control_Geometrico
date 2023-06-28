function Q = maxinv(A,X)
% MAXINV Q = maxinv(A,X) is an orthonormal basis for the maximal
%	 A-invariant contained in imX .
%	 The routine implements Relation 3.2.7 of "Controlled and
%	 Conditioned Invariant Subspaces in Linear System Theory"

%	 Basile and Marro 4-20-90

   Q = ortc(mininv(A',ortc(X)));
% --- last line of maxinv ---
