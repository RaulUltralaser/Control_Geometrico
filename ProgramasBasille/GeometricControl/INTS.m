function Q = ints(A,B)
% INTS  Intersection of subspaces Q = ints(A,B) is an orthonormal basis
%       for the subspace (imA) intersection (imB) .

%       Basile and Marro 4-20-90

   Q = ORTC(SUMS(ORTC(A),ORTC(B)));
% --- last line of ints
