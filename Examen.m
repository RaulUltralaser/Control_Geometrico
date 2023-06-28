clc
clearvars
% syms a b c d e f x y z c11 c12 c13 c21 c22 c23
% 
% A=[1 1 0; 0 2 0;0 0 0];
% C=[1 -1 0; 0 0 1];
% 
% O = obsv(A,C);
% 
% N = null(O);
% 
% Ab=[a b c;d e f; x y z];
% 
% Cb=[c11 c12 c13; c21 c22 c23];
% 
% P=[0 0 0; 0 0 0;0 0 1];

%% Sistema no observable 

A=[1 1;0 1];

B=[1;0];

C=[1 0];

O=obsv(A,C)

