syms s f1 f2
I=eye(2);
yd=5;

A=[-1 -2;
    0 -2];

B=[1;
   1];

S=[3;
   2];

D=[2 -3];

F1=[0 -1];

F2=[-3 4.5];

Fp=[-2/3 0];

P=[-2 -3];

K=place(A,B,P);

H=D*inv((s*I-A+B*K))*B;

H0=-.5;

Vd=yd/H0;


%% Sacar f
a=-4;
[S1, S2]=solve(3*f1+2*f2-4==(3/2)*a,3*f1+2*f2-2==a)

