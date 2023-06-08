syms x l

A=[0 1 0  0 0;
   0 0 0  0 0;
   0 0 -2 0 0;
   0 0 0  0 1;
   0 0 0  0 0];

B=[0 0;
   1 1;
   0 0;
   0 0;
   1 -1];
matriz=l*eye(5)-A;

smithForm(matriz)

[V,D]=eig(A);

P=poly2sym(charpoly(A),x)