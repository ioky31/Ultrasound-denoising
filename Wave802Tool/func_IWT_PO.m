function x = func_IWT_PO(wc,L,qmf)
% IWT_PO -- Inverse Wavelet Transform (periodized, orthogonal)
%  Usage
%    x = IWT_PO(wc,L,qmf)
%  Inputs
%    wc     1-d wavelet transform: length(wc) = 2^J.
%    L      Coarsest scale (2^(-L) = scale of V_0); L << J;
%    qmf    quadrature mirror filter
%  Outputs
%    x      1-d signal reconstructed from wc
%
%  Description
%    Suppose wc = FWT_PO(x,L,qmf) where qmf is an orthonormal quad. mirror
%    filter, e.g. one made by MakeONFilter. Then x can be reconstructed by
%      x = IWT_PO(wc,L,qmf)
%
%  See Also
%    FWT_PO, MakeONFilter
%
    wcoef = func_ShapeAsRow(wc);
	x = wcoef(1:2^L);
	[n,J] = func_dyadlength(wcoef);
	for j=L:J-1
		x = func_UpDyadLo(x,qmf) + func_UpDyadHi(wcoef(func_dyad(j)),qmf)  ;
	end
    x = func_ShapeLike(x,wc);

%
% Copyright (c) 1993. Iain M. Johnstone
%     
    
    
%   
% Part of WaveLab Version 802
% Built Sunday, October 3, 1999 8:52:27 AM
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail wavelab@stat.stanford.edu
%   
    
