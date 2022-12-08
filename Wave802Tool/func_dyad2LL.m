function [ix,iy] = func_dyad2LL(j)
%
%%%%%%%It's My Own Function!!!!
%
% dyad2LL-- Index entire j-th dyad of 2-d wavelet xform in LL
%  Usage
%    [ix,iy] = dyad2LL(j);
%  Inputs
%    j     integer
%  Outputs
%    [ix,iy]    list of all indices of wavelet coeffts at j-th level
%
    ix = 1:(2^(j)) ;
    iy = 1:(2^(j)) ;