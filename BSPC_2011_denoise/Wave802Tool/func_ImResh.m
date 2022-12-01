function [Im] = ImResh(X,J,I)
% Wavelet vector form to Wavelet image form transform
% [Im] = ImResh(X,J) returns the Wave-Image version
%  of the Wave-Vector :
%
%  X = [Wa(J) | Wh(J) | Wv(J) | Wd(J) | Wh(J-1) | Wv(J-1) | Wd(J-1) | ...
%       ... | Wh(1) | Wv(1) | Wd(1)]
%
%  Im = [
%         Wa(J) | Wh(J)  |             |           |                |
%         ---------------|   Wh(J-1)   |           |                |
%         Wv(J) | Wd(J)  |             |  Wh(J-2)  |                |
%         ---------------------------- |           | ...  Wh(1)     |
%            Wv(J-1)     |   Wd(J-1)   |           |                |
%         -----------------------------------------|                |
%                    Wv(J-2)           |  Wh(J-2)  |                |
%         ----------------------------------------------------------|
%                                .                 | .              |
%                                .                 |   .            |
%                                .                 |     .          |
%                                                  |                |
%                               Wv(1)              |      Wd(1)     |
%                                                  |                |
%         -----------------------------------------------------------
%        ]
%

[n,T] = size(X);
I0    = sqrt(T);
Im    = zeros(I0,I0);

% ---- Image lengths at each resolution

if nargin == 2
  I = DefineDyadLengths(I0,J);
end

% ---- First Resolution -------

len = I(J)^2;

Im(1:I(J),1:I(J))               = reshape(X(1:len),I(J),I(J));
Im(1:I(J),I(J)+1:2*I(J))        = reshape(X(len+1:2*len),I(J),I(J));
Im(I(J)+1:2*I(J),1:I(J))        = reshape(X(2*len+1:3*len),I(J),I(J));
Im(I(J)+1:2*I(J),I(J)+1:2*I(J)) = reshape(X(3*len+1:4*len),I(J),I(J));

% ---- Intermediate Resolution -------

for j = J-1:-1:1
  len = I(j)^2;
  
  Im(1:I(j),I(j)+1:2*I(j))        = reshape(X(len+1:2*len),I(j),I(j));
  Im(I(j)+1:2*I(j),1:I(j))        = reshape(X(2*len+1:3*len),I(j),I(j));
  Im(I(j)+1:2*I(j),I(j)+1:2*I(j)) = reshape(X(3*len+1:4*len),I(j),I(j));
end
