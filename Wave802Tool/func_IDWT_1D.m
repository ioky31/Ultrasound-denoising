function [X,scale] = func_IDWT_1D(W,resol,qmf)
% Inverse 1D-wavelet transform
% [X,scale] = IWaveletTransform(W,resolution,qmf)
%  1D-inverse wavelet transform: returns in X[n,T] the 
%  inverse 1D-wavelet transform of W[n,T] along the second 
%  dimension (i.e: each row represents a 1D signal).
%
% INPUTS
% ------
%  W          : [n by T] matrix of 1D-wavelet coefficients
%  resolution : resolution level,
%  qmf        : quadratic mirror filters.
%
% DEFAULT: resolution = log2(T), qmf <=> {Symmlet,4}.
%
% OUTPUTS
% -------
%  X     : [n by T] matrix of 1D signals
%  scale : scale parameter (ineteger)
%
% see also: WaveletTransform, MakeONFilter, IFWT_PO
%

scale = log2(size(X,2))-resol; 

X = func_IWT_PO(W,scale,qmf);



% [m,T] = size(W);
% 
% if nargin <= 2
%   qmf = func_MakeONFilter('Symmlet',4);
%   if nargin <= 1
%     resolution = log2(T);
%   end
% end
% 
% scale = log2(size(X,2))-resol; 
% 
% scale = log2(T) - resolution;
% 
% for i = 1:m
%   X(i,:) = func_IWT_PO(W(i,:),scale,qmf);
% end
