function [ix,iy] = func_dyad2LH(j)
%
%%%%%%%It's My Own Function!!!!
%
% dyad2LH-- Index entire j-th dyad of 2-d wavelet xform in LH(right-top corner)
%  Usage
%    [ix,iy] = dyad2LH(j);
%  Inputs
%    j     integer
%  Outputs
%    [ix,iy]    list of all indices of wavelet coeffts at j-th level
%
    ix = 1:(2^j) ;
    iy = (2^(j)+1):(2^(j+1)) ;