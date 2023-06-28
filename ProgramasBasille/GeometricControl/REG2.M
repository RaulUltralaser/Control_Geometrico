%   REG2 - Regulator Design
%
%   This program presents step by step the design of a standard autonomous
%   regulator with the Francis algorithm .
%   The program requires the following matrices, which must be defined
%   in an m-file whose name will be specified in interactive mode :
%
%   A1, B1, A2, A3, E1, E2, J .
%
%   The free eigenvalues of the plant and of the observer can be defined
%   in the same file or are entered in interactive mode .
%
%   At the end of the design, it is possible to perform a robustness
%   check, where matrices A1, B1, A3, L and M are varied at random by an
%   arbitrary percentage which is communicated in interactive mode .
%   Variation can affect either all the elements or only the nonzero
%   elements of the above matrices .
%   Graphic output is recorded in file reg2.met .
%   This program implements Algorithm 6.3-3 of "Controlled and Conditioned
%   Invariant Subspaces in Linear System Theory" .

%   Basile and Marro 4-20-90

clear
!del reg2.met
xx = '   strike any key to continue';
help reg2;
nprog = 2;
nopt = 0;
disp(xx)
pause
disp(' ')
disp('   enter the name of your data file; files certainly available are :')
disp('   examp1, examp2, examp3, examp4')
name = input('   name : ','s');
eval(name)
disp(xx)
pause
disp(' ')
int = 0;
  while (int ~= 1)&(int ~= 2)
  disp(' ')
  disp('   1 - use eigenvalues predefined in column vectors P1 and P2')
  disp('   2 - enter the free eigenvalues in interactive mode')
  int = input('enter your choice : ');
end
idis = 0;
  while (idis ~= 1)&(idis ~= 2)
  disp(' ')
  disp('   1 - no information on changes of basis')
  disp('   2 - information on changes of basis')
  idis = input('enter your choice : ');
end
idis = idis-1;
disp(' ')
disp('   matrix A1')
A1
disp(xx)
pause
disp(' ')
disp('   matrix A2')
A2
disp(xx)
pause
disp(' ')
disp('   matrix A3')
A3
disp(xx)
pause
disp(' ')
disp('   matrix B1')
B1
disp(xx)
pause
disp(' ')
disp('   matrix E1')
E1
disp(xx)
pause
disp(' ')
disp('   matrix E2')
E2
disp(xx)
pause
disp(' ')
disp('   matrix J')
J
disp(xx)
pause
disp(' ')
[t,na1] = size(A1);
[t,na2] = size(A2);
[me,t] = size(E1);
[t,nb] = size(B1);
[t,nJ] = size(J);
A4 = zeros(na2,na1);
B2 = zeros(na2,nb);
disp('   matrix A')
A=[A1 A3; A4 A2];
[ma,na] = size(A);
A
disp(xx)
pause
disp(' ')
disp('   matrix B')
B=[B1; B2];
B
disp(xx)
pause
disp(' ')
disp('   matrix E')
E=[E1 E2];
E
disp(xx)
pause
disp(' ')
disp('   poles of the plant')
eig(A1)
disp(xx)
pause
disp(' ')
disp('   invariant zeros of the plant : must be different from the')
disp('   eigenvalues of the exosystem to guarantee complementability')
gazero(A1,B1,E1)
disp(xx)
pause
disp(' ')
kerE=ker(E);
kerE
disp(xx)
pause
disp(' ')
disp('   maximal (A,B)-controlled invariant contained in kerE')
Vm=mainco(A,B,kerE);
Vm
disp(xx)
pause
disp(' ')
disp('   plant')
Pl=[eye(na1); zeros(na2,na1)];
Pl
disp(xx)
pause
disp(' ')
disp('   check if Vm complements the plant')
disp(' ')
disp('   Vm + Pl - must be the whole space')
sums(Vm,Pl)
disp(xx)
pause
disp(' ')
disp('   Vm ints Pl')
ints(Vm,Pl)
disp(xx)
pause
disp(' ')
disp('   construction of the regulator')
disp(' ')
disp('   matrix A1e')
A1e=A1;
A1e
disp(xx)
pause
disp(' ')
disp('   matrix A2e (replicated internal model)')
na2e = me*nJ;
A2e = zeros(na2e);
for i = 1:me
  i1=((i-1)*nJ+1);
  i2=i*nJ;
  A2e(i1:i2,i1:i2)=J(1:nJ,1:nJ);
end
A2e
disp(xx)
pause
disp(' ')
disp('   matrix A3e (arbitrary such that (Ae,Ee) is observable)')
disp(' ')
disp('   search for an observability path')
A3e = zeros(na1,na2e);
[t,I]=sort(abs(E1'));
for i = 1:me
  i1=((i-1)*nJ+1);
  i2=I(na1,i);
  A3e(i2,i1)=1;
end
A3e
disp(xx)
pause
disp(' ')
disp('   matrix A4e')
A4e = zeros(na2e,na1);
A4e
disp(xx)
pause
disp(' ')
disp('   matrix Ae')
Ae = [A1e A3e; A4e A2e];
[t,nae]=size(Ae);
Ae
disp(xx)
pause
disp(' ')
disp('   matrix Be')
Be = [B1; zeros(na2e,nb)];
Be
disp(xx)
pause
disp(' ')
disp('   matrix Ee')
Ee=[E1 zeros(me,na2e)];
Ee
disp(xx)
pause
disp(' ')
disp('   maximal A2e-invariant contained in kerA3e (must be the origin)')
maxinv(A2e,ker(A3e))
disp(xx)
pause
disp(' ')
kerEe=ker(Ee);
disp('   maximal Ae-invariant contained in kerEe (must be the origin)')
maxinv(Ae,kerEe)
disp(xx)
pause
disp(' ')
disp('   the plant in the regulator state space')
Ple = [eye(na1); zeros(na2e,na1)];
Ple
disp(xx)
pause
disp(' ')
disp('   maximal (Ae,Be)-controlled invariant contained in kerEe')
Vme=mainco(Ae,Be,kerEe);
Vme
disp(xx)
pause
disp(' ')
disp('   check if Vme complements Ple')
disp(' ')
disp('   Vme + Ple -must be the whole space')
sums(Vme,Ple)
disp(xx)
pause
disp(' ')
disp('   Vme ints Ple')
ints(Vme,Ple)
disp(xx)
pause
disp(' ')
tol=norm(A,'fro')*eps*10^4;
Rp=ints(Vme,Ple);
[t,nvi]=size(Rp);
no=norm(Rp,'fro');
if (nvi == 1)&(no < tol)
  nvi = 0;
end
if nvi ~= 0
  disp('   since Vme is greater than a complement of the plant')
  disp('   we construct from it a complement of the plant')
  disp(' ')
  Ve = compl(Ae,Be,Ee,Ple,idis);
else
  Ve = Vme;
end
disp('   the resolvent')
  Ve
  disp(xx)
  pause
  disp(' ')
disp('   checks about Ve')
disp(' ')
disp('   Ve + Ple - must be the whole space')
sums(Ve,Ple)
disp(xx)
pause
disp(' ')
disp('   Ve ints Ple - must be the origin')
ints(Ve,Ple)
disp(xx)
pause
disp(' ')
disp('   Ee * Ve - must be zero')
Ee*Ve
disp(xx)
pause
disp(' ')
disp('   pole assignment for the plant and derivation of an Fe such')
disp('   that (Ae+Be*Fe)*Ve is contained in Ve')
disp(' ')
if int == 2
  Fe = stabf(Ae,Be,Ve,Ple,idis);
else
  disp('   eigenvalues to be assigned')
  P1t=P1(1:na1,:);
  P1t
  disp(xx)
  pause
  disp(' ')
  Fe = stabf(Ae,Be,Ve,Ple,idis,P1t);
end
disp(' ')
disp('   matrix Fe')
Fe
disp(xx)
pause
disp(' ')
disp('   pole assigment for the regulator')
disp(' ')
if int == 2
  Ge = place1(Ae',Ee');
  Ge = Ge';
else
  disp('   eigenvalues to be assigned')
  P2t=P2(1:nae,:);
  P2t
  disp(xx)
  pause
  disp(' ')
  Ge = place1(Ae',Ee',P2t);
  disp(' ')
  Ge = Ge';
end
disp(' ')
disp('   matrix Ge')
Ge
disp(xx)
pause
disp(' ')
disp('   construction of the regulator')
disp(' ')
disp('   matrix A1hat = A')
A1hat=A;
A1hat
disp(xx)
pause
disp(' ')
disp('   matrix A3hat = B*Fe')
A3hat=B*Fe;
A3hat
disp(xx)
pause
disp(' ')
disp('   matrix A4hat = -Ge*E')
A4hat=-Ge*E;
A4hat
disp(xx)
pause
disp(' ')
disp('   matrix A2hat = Ae+Be*Fe+Ge*Ee')
A2hat=Ae+Be*Fe+Ge*Ee;
A2hat
disp(xx)
pause
disp(' ')
disp('   eigenvalues of A2hat')
eig(A2hat)
disp(xx)
pause
disp(' ')
%disp('   matrix Ahat')
Ahat=[A1hat A3hat; A4hat A2hat];
%ahat
%disp(xx)
%pause
%disp(' ')
%disp('   matrix Ehat')
Ehat=[E zeros(me,nae)];
%Ehat
%disp(xx)
%pause
%disp(' ')
disp('   matrix Aplant (matrix A of the regulated plant)')
A1plant=Ahat(1:na1,1:na1);
A3plant=Ahat(1:na1,(na+1):(na+nae));
A4plant=Ahat((na+1):(na+nae),1:na1);
A2plant=Ahat((na+1):(na+nae),(na+1):(na+nae));
Aplant=[A1plant A3plant; A4plant A2plant];
Aplant
disp(xx)
pause
disp(' ')
disp('   eigenvalues of Aplant (must be stable)')
eig(Aplant)
disp(xx)
pause
disp(' ')
disp('   matrix Bplant (matrix B of the regulated plant)')
Bplant=[Ahat(1:na1,(na1+1):na) ; Ahat((na+1):(na+nae),(na1+1):na)];
Bplant
disp(xx)
pause
disp(' ')
disp('   matrix Cplant (matrix C of the regulated plant)')
Cplant=[E1 zeros(me,nae)];
Cplant
disp(xx)
pause
disp(' ')
disp('   matrix Dplant (matrix D of the regulated plant)')
Dplant=E2;
Dplant
disp(xx)
pause
disp(' ')
disp('   computation of the invariant via the Sylvester equation')
%A=aplant;
%B=ahat(4:6,4:6);
%C=[ahat(1:3,4:6) ; ahat(7:12,4:6)];
%xin=lyap(A,-B,C);
Xin=lyap(Aplant,-Ahat((na1+1):na,(na1+1):na),Bplant);
X1=Xin(1:na1,:);
X2=Xin((na1+1):(na1+na1),:);
X3=Xin((na1+na1+1):(na1+nae),:);
Wsyl=[X1;eye(na2);X2;X3];
Wsyl
disp(xx)
pause
disp(' ')
disp('   Ehat*Wsyl - must be zero')
Ehat*Wsyl
disp(xx)
pause
disp(' ')
disp('   internal eigenvalues of Wsyl - must include those of the exosystem')
[P,Q]=stabi(Ahat,Wsyl);
eig(P)
disp(xx)
pause
disp(' ')
disp('   external eigenvalues of Wsyl - must be all stable')
eig(Q)
disp(xx)
pause
disp(' ')
%disp('   invariant zeros')
%gazero(Aplant,Bplant,Cplant,Dplant)
%disp(xx)
%pause
%disp(' ')
disp('   maximal Ahat-invariant contained in kerEhat')
kerEhat=ker(Ehat);
W=maxinv(Ahat,kerEhat);
W
disp(xx)
pause
disp(' ')
disp('   internal eigenvalues of W - must include those of the exosystem')
[P,Q]=stabi(Ahat,W);
eig(P)
disp(xx)
pause
disp(' ')
disp('   external eigenvalues of W - must be stable')
eig(Q)
disp(xx)
pause
disp(' ')
disp('   check of congruence')
disp(' ')
disp('   Wsyl + W (must have the same dimension as Wsyl and W)')
sums(Wsyl,W)
disp(xx)
pause
disp(' ')
disp('   Wsyl ints W (must have the same dimension as Wsylv and W)')
ints(Wsyl,W)
disp(xx)
pause
disp(' ')
nrloc = 1;
axis(TA);
disp('  multivariable root locus')
Aloc = [A1 B1*Fe; zeros(nae,na1) A2hat];
Bloc = [zeros(na1,me); Ge];
Cloc = [E1 zeros(me,nae)];
Dloc = zeros(me);
disp(' ')
disp('   open-loop poles')
olpo = eig(Aloc);
olpo
disp(xx)
pause
disp(' ')
%disp('   open-loop zeros')
%olze = gazero(Aloc,Bloc,Cloc);
%olze
disp('   open-loop zeros (also those at infinity)')
kz=10^10;
olze = mlocus(Aloc,Bloc,Cloc,Dloc,kz);
olze=olze.';
format short e
olze
format
disp(xx)
pause
disp(' ')
disp('   waiting for root locus computation')
disp(' ')
k=logspace(-2,2,100);
r = mlocus(Aloc,Bloc,Cloc,Dloc,k);
r1 = mlocus(Aloc,Bloc,Cloc,Dloc,1);
while nrloc == 1
  disp('   symbols :')
  disp('   x = open-loop poles')
  disp('   o = open-loop zeros')
  disp('   + = nominal poles')
disp(' ')
disp(xx)
pause
  x1 = real(olpo);
  y1 = imag(olpo);
  x2 = real(olze);
  y2 = imag(olze);
  plot(real(r),imag(r),'.',real(r1),imag(r1),'+',x1,y1,'x',x2,y2,'o'), pause
  clc
nrloc=3;
  disp('do you want to see the root locus with different axis scaling ?')
  while (nrloc ~= 0)&(nrloc ~= 1)
    nrloc=input('enter 1 for yes, 0 for no : ');
  end
  if nrloc == 1
    TA=input('type [xmin xmax ymin ymax] : ');
    axis(TA);
    disp(' ')
  end
end
meta reg2
disp(' ')
disp('   Do you want to check robustness ? 1 for yes, 0 for not')
nit = 2;
while (nit ~= 0)&(nit ~= 1)
  nit = input('enter your choice : ');
end
if nit == 0
  disp(' ')
  return
end
%nit = 1;
while nit == 1
disp(' ')
disp('   robustness check')
disp(' ')
disp('enter coefficient for robustness check (maximal variation in per cent)')
coe=input('   coefficient : ');
coe=coe/100;
disp(' ')
ind = 0;
  while (ind ~= 1)&(ind ~= 2)
  disp('   1 - Vary all the elements of the matrices')
  disp('   2 - Vary only the nonzero elements of the matrices')
  ind = input('enter your choice : ');
end
disp(' ')
disp('   matrix A1r')
A1r = randma(A1,coe,ind);
A1r
disp(xx)
pause
disp(' ')
disp('   matrix A3r')
A3r = randma(A3,coe,ind);
A3r
disp(xx)
pause
disp(' ')
disp('   matrix B1r')
B1r = randma(B1,coe,ind);
B1r
disp(xx)
pause
disp(' ')
disp('   matrix Fer')
Fer = randma(Fe,coe,ind);
Fer
disp(xx)
pause
disp(' ')
disp('   matrix Ger')
Ger = randma(Ge,coe,ind);
Ger
disp(xx)
pause
disp(' ')
disp('   matrix Ar')
Ar=[A1r A3r; A4 A2];
Ar
disp(xx)
pause
disp(' ')
disp('   matrix Br')
Br=[B1r; B2];
Br
disp(xx)
pause
disp(' ')
disp('   construction of the regulator')
disp(' ')
disp('   matrix A1hat = Ar')
A1hat=Ar;
A1hat
disp(xx)
pause
disp(' ')
disp('   matrix A3hat = Br*Fer')
A3hat=Br*Fer;
A3hat
disp(xx)
pause
disp(' ')
disp('   matrix A4hat = -Ger*E')
A4hat=-Ger*E;
A4hat
disp(xx)
pause
disp(' ')
disp('   matrix A2hat = Ae+Be*Fe+Ge*Ee')
A2hat=Ae+Be*Fe+Ge*Ee;
A2hat
disp(xx)
pause
disp(' ')
disp('   eigenvalues of A2hat')
eig(A2hat)
disp(xx)
pause
disp(' ')
%disp('   matrix Ahat')
Ahat=[A1hat A3hat; A4hat A2hat];
%ahat
%disp(xx)
%pause
%disp(' ')
%disp('   matrix Ehat')
Ehat=[E zeros(me,nae)];
%Ehat
%disp(xx)
%pause
%disp(' ')
disp('   matrix Aplant (matrix A of the regulated plant)')
A1plant=Ahat(1:na1,1:na1);
A3plant=Ahat(1:na1,(na+1):(na+nae));
A4plant=Ahat((na+1):(na+nae),1:na1);
A2plant=Ahat((na+1):(na+nae),(na+1):(na+nae));
Aplant=[A1plant A3plant; A4plant A2plant];
Aplant
disp(xx)
pause
disp(' ')
disp('   eigenvalues of Aplant (must be stable)')
eig(Aplant)
disp(xx)
pause
disp(' ')
disp('   matrix Bplant (matrix B of the regulated plant)')
Bplant=[Ahat(1:na1,(na1+1):na) ; Ahat((na+1):(na+nae),(na1+1):na)];
Bplant
disp(xx)
pause
disp(' ')
disp('   matrix Cplant (matrix C of the regulated plant)')
Cplant=[E1 zeros(me,nae)];
Cplant
disp(xx)
pause
disp(' ')
disp('   matrix Dplant (matrix D of the regulated plant)')
Dplant=E2;
Dplant
disp(xx)
pause
disp(' ')
disp('   computation of the invariant via the Sylvester equation')
Xin=lyap(Aplant,-Ahat((na1+1):na,(na1+1):na),Bplant);
X1=Xin(1:na1,:);
X2=Xin((na1+1):(na1+na1),:);
X3=Xin((na1+na1+1):(na1+nae),:);
Wsylr=[X1;eye(na2);X2;X3];
Wsylr
disp(xx)
pause
disp(' ')
disp('   Ehat*Wsylr - must be zero')
Ehat*Wsylr
disp(xx)
pause
disp(' ')
disp('   internal eigenvalues of Wsylr - must include those of the exosystem')
[P,Q]=stabi(Ahat,Wsylr);
eig(P)
disp(xx)
pause
disp(' ')
disp('   external eigenvalues of Wsylr - must be all stable')
eig(Q)
disp(xx)
pause
disp(' ')
%disp('   invariant zeros')
%gazero(Aplant,Bplant,Cplant,Dplant)
%disp(xx)
%pause
%disp(' ')
disp('   maximal Ahat-invariant contained in kerEhat')
Wr=maxinv(Ahat,kerEhat);
Wr
disp(xx)
pause
disp(' ')
disp('   internal eigenvalues of Wr - must include those of the exosystem')
[P,Q]=stabi(Ahat,Wr);
eig(P)
disp(xx)
pause
disp(' ')
disp('   external eigenvalues of Wr - must be stable')
eig(Q)
disp(xx)
pause
disp(' ')
disp('   check of congruence')
disp(' ')
disp('   Wsylr + Wr (must have the same dimension as Wsylr and Wr)')
sums(Wsylr,Wr)
disp(xx)
pause
disp(' ')
disp('   Wsylr ints Wr (must have the same dimension as Wsylr and Wr)')
ints(Wsylr,Wr)
disp(xx)
pause
disp(' ')
%
disp('   do you want to check robustness of the pole layout ?')
nlay = 3;
  while (nlay ~= 0)&(nlay ~= 1)
    nlay=input('enter 1 for yes, 0 for no : ');
  end
if nlay == 1
  axis(TA);
  disp(' ')
  disp('   initial layout :')
  disp('   x = open-loop poles')
  disp('   o = open-loop zeros')
  disp('   + = nominal poles')
disp(' ')
disp(xx)
pause
  plot(real(r1),imag(r1),'+',x1,y1,'x',x2,y2,'o'), pause
  meta reg2
  nit = 1;
  while nit == 1;
    np = input('enter the number of points to compute : ');
    disp(' ')
    disp('   thinking')
    disp(' ')
    r = sqrt(-1)*ones((na1+nae),np);
    A2er = Ae+Be*Fe+Ge*Ee;
    for i = 1:np
      A1r=randma(A1,coe,ind);
      B1r=randma(B1,coe,ind);
      Fer=randma(Fe,coe,ind);
      Ger=randma(Ge,coe,ind);
      Ar = [A1r B1r*Fer;-Ger*E1 A2er];
      r(:,i) = eig(Ar);
    end
    r = r.';
    plot(real(r),imag(r),'.',real(r1),imag(r1),'+',x1,y1,'x',x2,y2,'o'), pause
    clc
    nit=3;
    disp('   do you want more points ?')
    while (nit ~= 0)&(nit ~= 1)
      nit=input('enter 1 for yes, 0 for no : ');
    end
  hold on
  end
  hold off
  meta reg2
end
nit=3;
disp(' ')
disp('do you want to repeat the robustness check ?')
  while (nit ~= 0)&(nit ~= 1)
    nit=input('enter 1 for yes, 0 for no : ');
  end
disp(' ')
end
% --- last line of reg2 ---
