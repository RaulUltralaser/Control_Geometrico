clc
clearvars

%%%% sistema
A=[2 1 0;0 -3 0;-1 0 -1];
B=[0;1;0];
C=[0 1 0; 0 0 1];
L1=[0;1;2];
L=[0 ;0; 1];  %tecnicamente este el L2, pero para el observador lo defini 
              %como L, así que prefiero dejarlo así para que siga
              %funcionando


%% Esto es lo para el cálculo de w_estrella
kerC=null(C,'r');     %inicializó el kernel de C
w_prev=L1;             %inicializó e

%% Este ciclo es la iteración para w_estrella
for i= 1:10
    inters=intersection(w_prev,kerC); %usa la fución intersection para calular la intersection de w_i-1 con KerC
    w_i=cat(2,w_prev,A*inters)       %Los suma como nos dijo Ruíz no es suma directamente sino concatenación
    w_prev                           %Esto lo dejo para debugguear más fácil (que me imprima el valor anterior)
%     if isequal(w_i(:,i+1),w_i(:,i))   %debería revisar que la última columna sea linealmente independiente no iguales
%         w_i=w_i(:,1:i);               %pero por ahora así jala
%     else
%         w_i=w_i;
%     end
    
    if isequal(w_i,w_prev) %compara solo la última columna de cada w, ya que se concatenan
        disp('w_estrella=')
        w_i                       %si son iguales lo imprime
        break
    else
        w_prev=w_i;                      %si No son iguales ahora w_i es w_i-1
    end
    i=i+1;
end 
%% Esto es para calcular la representación matricial de la proyección canonica
syms b11 b12 b21 b22 b31 b32
P=null(w_i','r')'; %P*W_i deben ser sero, por lo que así está bien

D=[b11 b12;
   b21 b22;
   b31 b32];    %Defino D simbolicamente

MP=P*(A+D*C)*w_i; % De esta ecuación obtenemos que mientras b22 y b12 sea igual a cero está bien

D=[1 0;
   1 0;
   1 1];       %Defino D numericamente

P*(A+D*C)*w_i; %Esto es para saber si está ecuación es igual a cero


Pr=[1 0;
    0 1;
    0 0];      %P*Pr debe ser igual a I
P*Pr;          %Esto es para confirmar que sea igual a I
%% Defino cosas del observador
%Zdot=P*A0*Pr+P*B*u-P*D*y
%g=C*Pr
A0=A+D*C;      

F=P*A0*Pr;      
G=P*B;
E=P*D;
M=C*Pr;

%% Esto es para una ley de control que vuelva estable mi sistema
R=place(A,B,[-1 -2 -3]);

%% Esto es lo para el cálculo de S_estrella
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

%% Defino cosas para el generador de residuos
%H es el aniquilador de CS (NOTAS DEL 26 de JULIO)
% H=null(C*S_i,'r');      %Valio verga porque no tiene kernel
syms H11 H12 H21 H22
Hs=[H11 H12;H21 H22];
% Hice los cálculos  a mano en mi libretita
H=[2 -1;2 -1];


%% Función para calcular la intersección de subespacios
function Cap = intersection(A,B)
    Aort = null(A','r')';
    Bort = null(B','r')';
    Suma = cat(1, Aort, Bort);
    Cap = null(Suma,'r');   
end



