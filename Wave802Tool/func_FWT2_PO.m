function wc = func_FWT2_PO(x,L,qmf)
% FWT2_PO -- 2-d MRA wavelet transform (periodized, orthogonal)
%  Usage
%    wc = FWT2_PO(x,L,qmf)
%  Inputs
%    x     2-d image (n by n array, n dyadic)
%    L     coarse level
%    qmf   quadrature mirror filter
%  Outputs
%    wc    2-d wavelet transform
%
%  Description
%    A two-dimensional Wavelet Transform is computed for the
%    array x.  To reconstruct, use IWT2_PO.
%
%  See Also
%    IWT2_PO, MakeONFilter
%
	[n,J] = func_quadlength(x);
	wc = x; 
	nc = n;
	for jscal=J-1:-1:L,
		top = (nc/2+1):nc; bot = 1:(nc/2);
		for ix=1:nc,
			row = wc(ix,1:nc);
			wc(ix,bot) = func_DownDyadLo(row,qmf);
			wc(ix,top) = func_DownDyadHi(row,qmf);
		end
		for iy=1:nc,
			row = wc(1:nc,iy)';
			wc(top,iy) = func_DownDyadHi(row,qmf)';
			wc(bot,iy) = func_DownDyadLo(row,qmf)'; 
		 end
		nc = nc/2;
	end   
 
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
    
