%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%NOTA: Es necesario correr w_estrella primero%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars -except w_i

%%%% sistema
A=[2 1 0;0 -3 0;-1 0 -1];
B=[0;1;0];
C=[0 1 0; 0 0 1];
L=[0 ;0; 1];

%% Esto es lo para el cálculo de S_estrella
kerC=null(C,'r');     %inicializó el kernel de C
X=eye(size(A));      %inicializó el espacio X base canónica
S_prev=X;             %inicializó s0


while i<10
    P1=inv(A)*S_prev;
    P2=intersection(P1,kerC);
    S_i=cat(2,w_i,P2);
    if isequal(S_i,S_prev)
        disp('S_estrella es:')
        S_i
        break
    else
        S_prev=S_i;
    end
    i=i+1;
end 



function Cap = intersection(A,B)
    Aort = null(A','r')';
    Bort = null(B','r')';
    Suma = cat(1, Aort, Bort);
    Cap = null(Suma,'r');   
end