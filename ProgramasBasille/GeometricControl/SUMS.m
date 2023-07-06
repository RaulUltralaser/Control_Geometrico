function Q = sums(A,B)
% SUMS  Sum of subspaces Q = sums(A,B) is an orthonormal basis for
%       the subspace im[A B] = imA + imB .

%       Basile and Marro 4-20-90

   Q = IMA([A B],1);
% --- last line of sums ---
