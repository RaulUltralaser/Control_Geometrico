function B = randma(A,coe,ik)
% RANDMA B = randma(A,coe,ik)
%	 Random variation of a matrix. It is used in programs REG* to
%	 perform robustness checks. If ik = 1 all elements of matrix
%	 A are subject to vary with respect to their nominal values
%	 in an interval with length coe, if ik = 2 only those which
%	 are different from zero are subject to vary by a percentage
%	 of their values.

%	 Basile and Marro 4-20-90

[ma,na] = size(A);
if (ik ~= 1)&(ik ~= 2)
  error('input index wrong in RANDMA');
end
if ik == 2
  B=A;
  for i = 1:ma
    for j = 1:na
      B(i,j) = A(i,j)+coe*A(i,j)*(1-2*rand(1));
    end
  end
else
  B=A+coe*(ones(ma,na)-2*rand(ma,na));
end
% --- last line of randma ---
