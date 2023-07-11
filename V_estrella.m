clc
clearvars

%% Se tienen que agregar las matrices manualmente

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

D=[0 0 0 1 0];
%% Esto asigna los valores que necesito antes de empezar con el ciclo
KerD=null(D,'r');
V_i = KerD;
V_prev = V_i;
r = size(KerD);
%% EL ciclo For se detiene cuando un V_i es igual a V_{i+1}, y lo muestra
%%en la command window
for i= 0:r
    w_i=null(cat(2,B,V_prev)','r')';
    cat(1,D,w_i*A);

    V_i=null(cat(1,D,w_i*A),'r');
    if isequal(V_i,V_prev)
        V_i
        break
    else
        V_prev=V_i;
    end
    i=i+1;
end




%% Esta función no la uso aquí pero podría ser util, calcula la intersección de dos conjuntos
function Cap = intersection(A,B)
    Aort = null(A','r')';
    Bort = null(B','r')';
    Suma = cat(1, Aort, Bort);
    Cap = null(Suma,'r');   
end

