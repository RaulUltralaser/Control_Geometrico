function Q = invt(A,X)
% INVT	 Inverse transform. Q=invt(A,p) is an orthonormal basis for the
%	 inverse map of imX in A.

%	 Basile and Marro 4-20-90

Q=ortc(A'*ortc(X))
% --- last line of invt ---
