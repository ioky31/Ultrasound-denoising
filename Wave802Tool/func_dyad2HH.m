function [ix,iy] = func_dyad2HH(j)
%
%%%%%%%It's My Own Function!!!!
%
% dyad2HH-- Index entire j-th dyad of 2-d wavelet xform in HH
%  Usage
%    [ix,iy] = dyad2HH(j);
%  Inputs
%    j     integer
%  Outputs
%    [ix,iy]    list of all indices of wavelet coeffts at j-th level
%
    ix = (2^(j)+1):(2^(j+1)) ;
    iy = (2^(j)+1):(2^(j+1)) ;