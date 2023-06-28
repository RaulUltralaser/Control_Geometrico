%       EXAMP1 is a file of data for REG1, REG2 or REG3
%       The data are matrices A1, A2, A3, B1, E1 and E2. Poles to be
%       assigned are specified in the column vector P1 for the plant and
%       P2 for the regulator, but can be also communicated in interactive
%       mode if preferred . A further datum for REG2 and REG3 are the
%       Jordan block J to be replicated in the exosystem, and for REG3
%       specifically, the zeros to be assigned and the loop gain.

% Matrices for REG1, REG2 and REG3
%
help examp1;
%
  A1=[-.4 -4.04 0; 1 0 0; 0 0 -1];
  B1=[2 0; .4 1; 0 1.5];
  E1=[0 -1 0; 0 0 -1];
  E2=[1 0; 0 1];
  A3=[0 0;.2 0;0 .3];
  A2=[0 1; 0 0];
%
% Different options for the exosystem (considered in REG1)
%
if nopt == 2
  A3=[0 0 0; .2 0 0; 0 0 .3];
  A2=[0 1 0; 0 0 0; 0 0 0];
  E2=[1 0 0; 0 0 1];
end
if nopt == 3
  A3=[0 0 0 0; .2 0 0 0; 0 0 .3 0];
  A2=[0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
  E2=[1 0 0 0; 0 0 1 0];
end
%
% Poles to be assigned (they can be changed in interactive mode)
%
  j=sqrt(-1);
  P1=[-1-j;-1+j;-2];  % poles for the plant
  P2=[-1-2*j;-1+2*j;-.8;-1.5;-3;-6;-7]; % poles for the regulator
%
%  Internal model structure to be repeated (for REG2 and REG3)
%
if (nprog == 2)|(nprog == 3)
  J=[0 1;0 0];
end
%
%  Zeros and initial value of the loop gain (for REG3)
%
if nprog == 3
  Z1=[-2-j;-2+j;-2+1.8*j;-2-1.8*j];
  gain = 20;
end
%
%  Set axis graduation for the root locus
%
TA = [-8 2 -10 10];
% --- last line of examp1 ---
