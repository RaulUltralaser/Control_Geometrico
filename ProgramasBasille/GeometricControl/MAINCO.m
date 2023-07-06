function Q = mainco(A,B,X)
% MAINCO Q = mainco(A,B,X) is an orthonormal basis for the maximal
%	 (A,imB)-controlled invariant contained in imX .
%	 The routine implements Relation 4.1.7 of "Controlled and
%	 Conditioned Invariant Subspaces in Linear System Theory"

%	 Basile and Marro 4-20-90

   Q = ORTC(MIINCO(A',ORTC(B),ORTC(X)));
% --- last line of mainco ---
