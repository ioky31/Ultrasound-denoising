function [W,resol,qmf] = func_DWT_1D(X,scale,wbase,mom)
% Forward 1D-wavelet transform
% [W,resolution,qmf] = WaveletTransform1D(X,scale,wbase,mom)
%  1D-wavelet transform: returns in W[n,T] the 1D-wavelet transform 
%  of X[n,T] along the second dimension (each row represents a 1D signal).
%
% INPUTS
% ------
%  X     : [n by T] matrix of 1D signals
%  wbase : wavelet basis ('Haar', 'Symmlet', 'Daubechies', ...)
%  scale : scale paramter (integer)
%  mom   : number of moments
%
% DEFAULT: wbase = 'Symmlet', mom = 4, scale = 0.
%
% OUTPUTS
% -------
%  W          : [n by T] matrix of 1D-wavelet coefficients
%  resolution : resolution level,
%  qmf        : quadratic mirror filters.
%
% see also: IWaveletTransform, MakeONFilter, FWT_PO
%

qmf = func_MakeONFilter(wbase,mom);

resol = log2(size(X,2))-scale; 

W = func_FWT_PO(X, scale, qmf);

% 
% [m,T] = size(X);
% W     = zeros(m,T);
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
% resolution = log2(T) - scale;
% 
% for i = 1:m
%   W(i,:) = func_FWT_PO(X(i,:),scale,qmf);
% end
