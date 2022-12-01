function [Xsub] = ReadSubMatrix(X,I1,subband)
% Wavelet subband extraction
% Xsub = ReadSubMatrix(X,I,subband returns in 
%  Xsub[n,I^2] the sub matrix of X[I,I,n] 
%  corresponding to a subband ('a', 'h', 'v' or 'd').
%

I2 = 2*I1;
switch subband
 case 'a'
  D = X(1:I1,1:I1,:);
 case 'h'
  D = X(1:I1,(I1+1):I2,:);
 case 'v'
  D = X((I1+1):I2,1:I1,:);
 case 'd'
  D = X((I1+1):I2,(I1+1):I2,:);
end

Xsub = D(:)';
