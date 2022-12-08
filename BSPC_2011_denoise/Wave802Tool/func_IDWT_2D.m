function [X,scale] = func_IDWT_2D(W,resol,qmf)
% Inverse 2D-wavelet transform
% [X,scale] = IWaveletTransform2D(W,resolution,qmf) inverse 
%  2D-wavelet transform: returns in X[n,T] the inverse
%  2D-wavelet transform of W[n,T] along the second dimension.
%  (i.e each row respresents an image)
%
% INPUTS
% ------
%  W          : [n by T] matrix of 2D-wavelet coefficients
%  resolution : resolution level,
%  qmf        : quadratic mirror filters.
%
% DEFAULT: resolution = log2(sqrt(T)), qmf <=> {Symmlet,4}.
%
% OUTPUTS
% -------
%  X     : [n by T] matrix of images
%  scale : scale parameter (ineteger)
%
% see also: WaveletTransform2D, MakeONFilter, IFWT2_PO
%

scale = log2(size(W,2)) - resol; 

Image  = func_IWT2_PO(W,scale,qmf);
X = Image(:)';


% [m,T] = size(W);
% I     = sqrt(T);
% Image = zeros(I,I);
% 
% if nargin <= 2
%   qmf = func_MakeONFilter('Symmlet',4);
%   if nargin <= 1
%     resolution = log2(I);
%   end
% end
% 
% scale = log2(I) - resolution;
% for i = 1:m
%   Image  = reshape(W(i,:),I,I);
%   Image  = func_IWT2_PO(Image,scale,qmf);
%   X(i,:) = Image(:)';
% end
