clc
clearvars

syms s a b c d e x y z w v e f g h i j k l m n o p
I=eye(5);

A=[0 0  0 0 0;
   0 -1 0 0 0;
   0 0  0 1 0;
   0 0  0 0 1;
   0 0  0 0 0];

B=[1 0;
   0 0;
   0 0;
   0 0;
   0 1];

S=[1 0;
   0 1;
   0 0;
   0 0;
   0 0];

D=[0 0 0 1 0];



P=[-1 -2 -3 -4 -5];
% 
% F=place(A,B,P);
F=[4 0 -13 -9 -1;
   0 0  30 31  10];

Af=A+B*F;

V_i=[1 0 0;0 1 0; 0 0 1; 0 0 0;0 0 0];

R_i=[1;0;0;0;];

P=[0 1 0 0 0;
   0 0 1 0 0;
   0 0 0 1 0;
   0 0 0 0 1];

Afbar=[a b c d;
       e f g h;
       i j k l;
       m n o p];
% Fsim=[a b c d e;x y z w v];
% 
% Abar=A+B*Fsim;


% obsv(A+B*Fsim,D)
% 
% Abar=A-B*F
% Bbar=[B zeros(5,3)]
% 
% [num, den] = ss2tf(Abar, Bbar, D, 0);

% H=D*inv((s*I-A+B*F))*B
% 
% H0=-.5;
% 
% Vd=yd/H0;