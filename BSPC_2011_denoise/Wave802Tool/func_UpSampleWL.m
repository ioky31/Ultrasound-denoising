function y = func_UpSample(x,s)
% UpSample -- Upsampling operator
%  Usage
%    u = UpSample(d[,s]) 
%  Inputs
%    d   1-d signal, of length n
%    s   upsampling scale, default = 2
%  Outputs
%    u   1-d signal, of length s*n with zeros
%        interpolating alternate samples
%        u(s*i-1) = d(i), i=1,...,n
%

	if nargin == 1, s = 2; end
	n = length(x)*s;
	y = zeros(1,n);
	y(1:s:(n-s+1) )=x;
    
    
%   
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
    
