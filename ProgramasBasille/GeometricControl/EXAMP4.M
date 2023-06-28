%       EXAMP4 is a file of data for REG1, REG2 or REG3
%       The data are matrices A1, A2, A3, B1, E1 and E2. Poles to be
%       assigned are specified in the column vector P1 for the plant and
%       P2 for the regulator, but can be also communicated in interactive
%       mode if preferred . A further datum for REG2 and REG3 are the
%       Jordan block J to be replicated in the exosystem, and for REG3
%       specifically, the zeros to be assigned and the loop gain.

% Matrices for REG1, REG2 and REG3
%
help examp4;
%
  A1=[-.4 -4.04 0; 1 0 0; 0 1 -1];
  B1=[1; 0; 0];
  E1=[0 0 -1];
  E2=[1 0];
%  E2=[1 0 0 0];
  A3=[0 0; 0 0; 0 .3];
%  A3=[0 0 0 0; 0 0 0 0; 0 .3 .2 .2];
  A2=[0 1; 0 0];
%  A2=[0 1 0 0; 0 0 0 0; 0 0 -6 0; 0 0 0 -5];
%
% Different options for the exosystem (considered in REG1)
%
if nopt == 2
  A3=[0 0 0; 0 0 0; 0 0 .3];
  A2=[0 1 0; 0 0 0; 0 0 0];
  E2=[1 0 0];
end
if nopt == 3
  A3=[0 0 0 0; 0 0 0 0; 0 0 .3 0];
  A2=[0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
  E2=[1 0 0 0];
end
%
% Poles to be assigned (they can be changed in interactive mode)
%
  j=sqrt(-1);
  P1=[-1-j;-1+j;-2];  % poles for the plant
  P2=[-1-2*j;-1+2*j;-1;-3;-5;-6;-7]; % poles for the regulator
%
%  Internal model structure to be repeated (for REG2 and REG3)
%
if (nprog == 2)|(nprog == 3)
%  J=[0 1;0 0];
  ina3 = 1;  %  to impose the following A3e
  A3e=[1 0 1 1; 0 0 0 0; 0 0 0 0];
  J=[0 1 0 0; 0 0 0 0; 0 0 -6 0; 0 0 0 -5];
end
%
%  Zeros and initial value of the loop gain (for REG3)
%
if nprog == 3
  Z1=[-.25-2*j;-.25+2*j;-.5+.3*j;-.5-.3*j];
  gain = 20;
end
%
%  Set axis graduation for the root locus
%
TA = [-8 2 -8 8];
% --- last line of examp4 ---
