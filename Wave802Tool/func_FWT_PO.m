function wcoef = func_FWT_PO(x,L,qmf)
% FWT_PO -- Forward Wavelet Transform (periodized, orthogonal)
%  Usage
%    wc = FWT_PO(x,L,qmf)
%  Inputs
%    x    1-d signal; length(x) = 2^J
%    L    Coarsest Level of V_0;  L << J
%    qmf  quadrature mirror filter (orthonormal)
%  Outputs
%    wc    1-d wavelet transform of x.
%
%  Description
%    1. qmf filter may be obtained from MakeONFilter   
%    2. usually, length(qmf) < 2^(L+1)
%    3. To reconstruct use IWT_PO
%
%  See Also
%    IWT_PO, MakeONFilter
%
  [n,J] = func_dyadlength(x) ;
  wcoef = zeros(1,n) ;
  beta = func_ShapeAsRow(x);  %take samples at finest scale as beta-coeffts
  for j=J-1:-1:L
       alfa = func_DownDyadHi(beta,qmf);
       wcoef(func_dyad(j)) = alfa;
       beta = func_DownDyadLo(beta,qmf) ;  
  end
  wcoef(1:(2^L)) = beta;
  wcoef = func_ShapeLike(wcoef,x);

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
    
