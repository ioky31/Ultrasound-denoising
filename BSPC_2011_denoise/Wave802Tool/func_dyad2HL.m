function [ix,iy] = func_dyad2HL(j)
%
%%%%%%%It's My Own Function!!!!
%
% dyad2HL-- Index entire j-th dyad of 2-d wavelet xform in HL(left-bottom corner)
%  Usage
%    [ix,iy] = dyad2HL(j);
%  Inputs
%    j     integer
%  Outputs
%    [ix,iy]    list of all indices of wavelet coeffts at j-th level
%
    ix = (2^(j)+1):(2^(j+1)) ;
    iy = 1:(2^j) ;