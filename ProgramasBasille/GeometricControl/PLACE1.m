function F=place1(A,B,P)
% PLACE1 Eigenvalue assignment
%        F = PLACE1(A,B,P) returns as F the state feedback matrix
%        such that the h assignable eigenvalues of A+B*F are those
%        specified in the first h rows of the column vector P .
%        F = PLACE1(A,B) assigns these eigenvalues in interactive mode .
%        To enter complex eigenvalues use symbol j for the imaginary unit .

%        Basile and Marro 4-20-90

j=sqrt(-1);
[ma,na]=size(A);
[mb,nb]=size(B);
T=mininv(A,B);
[m,n]=size(T);
nargs = nargin;
if n ~= 0
  if nargs == 2
  disp(   'Enter vector P of eigenvalues to be located: number of components')
    n
    ni=1;
    P=zeros(n,1);
    while (ni <= n)
      P(ni,1)=input(   'Enter an eigenvalue: ')
      ni=ni+1;
    end
  else
    if length(P) < n
      error('   not enough eigenvalues transmitted to PLACE1');
    end
    P=P(1:n,:);
  end
  disp('   thinking')
  disp(' ')
  if (n < na)
    T=ima([T ortc(T)],0);
    A1=inv(T)*A*T;
    B1=inv(T)*B;
    A1=A1(1:n,1:n);
    B1=B1(1:n,:);
    F1=place(A1,B1,P);
    F1=[F1 zeros(nb,(na-n))];
    F=-F1*T;
  else
    F=place(A,B,P);
    F=-F;
  end
else
  disp('   no eigenvalues can be assigned in PLACE1')
end
% --- last line of place1 ---
