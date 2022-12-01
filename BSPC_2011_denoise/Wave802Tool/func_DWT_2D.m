function [W,resol,qmf] = func_DWT_2D(X,scale,wbase,mom)
% Forward 2D-wavelet transform
% [W,resolution,qmf] = WaveletTransform2D(X,scale,wbase,mom)
%  2D-wavelet transform: returns in W[n,T] the 2D-wavelet transform
%  of X[n,T] along the second dimension (i.e each row represents an image).
%
% INPUTS
% ------
%  X     : [n by T] matrix of images
%  wbase : wavelet basis ('Haar', 'Symmlet', 'Daubechies', ...)
%  scale : scale parameter (ineteger)
%  mom   : number of vanishing moments
%
% DEFAULT: wbase = 'Symmlet', mom = 4, scale = 0.
%
% OUTPUTS
% -------
%  W          : [n by T] matrix of 2D-wavelet coefficients
%  resolution : resolution level,
%  qmf        : quadratic mirror filters.
%
% see also: IWaveletTransform2D, MakeONFilter, FWT2_PO
%

qmf = func_MakeONFilter(wbase,mom);

resol = log2(size(X,2)) - scale; 

Image  = func_FWT2_PO(X, scale, qmf);
W(i,:) = Image(:)';
  
% [m,T] = size(X);
% W     = zeros(m,T);
% I     = sqrt(T);
% Image = zeros(I,I);
% 
% %%% INPUT CHECK
% if nargin <= 3
%   mom = '';
%   if nargin <= 2
%     wbase = 'Symmlet';
%     mom   = 4;
%     if nargin <= 1
%       scale = 0;
%     end
%   end
% end
% 
% if isempty(mom)
%   qmf = func_MakeONFilter(wbase);
% else
%   qmf = func_MakeONFilter(wbase,mom);
% end
% %%% INPUT CHECK
% 
% resolution = log2(I) - scale;
% 
% for i = 1:m
%   Image  = reshape(X(i,:),I,I);
%   Image  = func_FWT2_PO(Image,scale,qmf);
%   W(i,:) = Image(:)';
% end
