function y = func_aconv(f,x)
% aconv -- Convolution Tool for Two-Scale Transform
%  Usage
%    y = aconv(f,x)
%  Inputs
%    f    filter
%    x    1-d signal
%  Outputs
%    y    filtered result
%
%  Description
%    Filtering by periodic convolution of x with the
%    time-reverse of f.
%
%  See Also
%    iconv, UpDyadHi, UpDyadLo, DownDyadHi, DownDyadLo
%

	n = length(x);
	p = length(f);
	if p < n,
	   xpadded = [x x(1:p)];
	else
	   z = zeros(1,p);
	   for i=1:p,
		   imod = 1 + rem(i-1,n);
		   z(i) = x(imod);
	   end
	   xpadded = [x z];
	end
	fflip = func_reverse(f);
	ypadded = filter(fflip,1,xpadded);
	y = ypadded(p:(n+p-1));

%
% Copyright (c) 1993. David L. Donoho
%     
    
    
%   
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
    
