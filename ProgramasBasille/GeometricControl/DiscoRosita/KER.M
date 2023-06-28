function Q = ker(A)
% KER   Kernel. Q=ker(A) is an orthonormal basis for kerA.

%       Basile and Marro 6-20-90

Q = ortc(A');
% --- last line of ker ---
