clc 
clearvars
close all 

%%%%%%%%%%%% Programa de minimo común divisor %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% calcula los divisores elementales de una matriz 3x3 %%%%%%%%%%%%%%
%%%% los presenta como delta=[delta1, delta2, delta3]%%%%%%%%%%%%%%%%%%%

%% Esto se tiene que poner manualmente 

A = [1 0 -1; 0 1 0; 0 0 2]; %%%Esta es la matriz a usar


%% Esto es para el cálculo de los divisores (no se pone nada manualmente)
syms l;                         %%% esto es lambda 
I=eye(3);                       %%% la matriz de unos

matriz = l*I-A;                 %%% lambda*I-A

%% esto es para calcular el primer divisor
%todo esta chingadera es porque gcd no me permite más de cuatro polinomios
% y ademas al usar 3 simbolicos tiene un error
d1=gcd(matriz(1,1),matriz(1,2));
d2=gcd(matriz(1,3),matriz(2,1));
d3=gcd(matriz(2,2),matriz(2,3));
d4=gcd(matriz(3,1),matriz(3,2));
dd1=gcd(d1,d2);
dd2=gcd(d3,d4);
ddd1=gcd(dd1,dd2);
delta1=gcd(ddd1,matriz(3,3));


%% esto es para calcular el segundo divisor
menores = calcularMenores(matriz);

d1=gcd(menores(1,1),menores(1,2));
d2=gcd(menores(1,3),menores(2,1));
d3=gcd(menores(2,2),menores(2,3));
d4=gcd(menores(3,1),menores(3,2));
dd1=gcd(d1,d2);
dd2=gcd(d3,d4);
ddd1=gcd(dd1,dd2);
delta2=gcd(ddd1,menores(3,3));

%% Esto es para calcular el tercer divisor

delta3=det(matriz);

%% presenta los divisores
delta = [delta1,delta2,delta3]

%% funciones que se mandan a llamar para el calculo del segundo divisor
% técnicamente separa la matriz en todas sus posibles combinaciones 
% y a cada combinación le calcula su determinante, acomoda los resultados
% en una matriz

function menores = calcularMenores(matriz)
    [n, ~] = size(matriz);
    menores = sym(zeros(n, n));
    
    for i = 1:n
        for j = 1:n
            menores(i, j) = calcularMenor(matriz, i, j);
        end
    end
end

function menor = calcularMenor(matriz, fila, columna)
    submatriz = matriz;
    submatriz(fila, :) = [];
    submatriz(:, columna) = [];
    menor = det(submatriz);
end

