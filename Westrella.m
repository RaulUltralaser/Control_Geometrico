clc
clearvars

%%%% sistema
A=[2 1 0;0 -3 0;-1 0 -1];
B=[0;1;0];
C=[0 1 0; 0 0 1];
e=[1 0 0];

kerC=null(C,'r');


w_prev=e;
% for i= 0:r2
%     wi=cat(2,A*S_prev,B);
%     S_i=intersection(V_i,wi);
%     if isequal(S_i,S_prev)
%         S_i
%         break
%     else
%         S_prev=S_i;
%     end
%     i=i+1;
% end 



function Cap = intersection(A,B)
    Aort = null(A','r')';
    Bort = null(B','r')';
    Suma = cat(1, Aort, Bort);
    Cap = null(Suma,'r');   
end
