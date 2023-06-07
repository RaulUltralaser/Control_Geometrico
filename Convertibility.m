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

[V,D]=eig(A);
syms x
P=poly2sym(charpoly(A),x)